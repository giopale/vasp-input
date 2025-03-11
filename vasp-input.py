import hydra
from pathlib import Path
from pymatgen.io.vasp import Poscar, Kpoints, Incar, Potcar, VaspInput
from omegaconf import OmegaConf
from omegaconf.listconfig import ListConfig
import numpy as np
import sys
from itertools import product
import json, yaml
import copy
import os, stat
import subprocess

from logger_setup import logger


def ask_if_overwrite(dd):
    folder = Path(dd)
    if folder.exists():
        user_input = input(f"{str(folder)} exists. Do you want to continue? (y/n): ")
        if user_input.lower() == "y":
            proceed = True
        else:
            proceed = False
    else:
        proceed = True
    return proceed


def decide_overwrite(dd, overwrite):
    if overwrite:
        return True
    else:
        return ask_if_overwrite(dd)


def is_single_line(aa):
    if isinstance(aa, str):
        lines = aa.strip().splitlines()
        if len(lines) == 1:  # Single-line string, likely a path
            return Path(aa).is_file()
    return False


def scan_vasp_files(keys, dict_files):
    # scan a directory for desired files
    # if specified, find other files and add them (with priority)
    # or manual input (higher priority)

    if dict_files.dir is None:
        logger.warning("source.dir not set")
    sourcedir = Path(dict_files.dir)
    results = {}

    # find files corresponding to keys
    for key in keys:
        value = next(sourcedir.glob(f"{key}"), None)
        results.update({key: value})

    # build out dict
    for key, value in dict_files.items():  # <--------------------------- to test
        if key == "dir":
            continue
        elif value is not None:
            if is_single_line(value):
                if Path(value).exists():
                    results.update({key.upper(): Path(value)})
                else:
                    logger.warning(f"{value} not found")
            else:
                results.update({key.upper(): value})
    err = False

    # if file provided manually, overwrite

    #    for key, val in results.items():
    #        if val is None:
    #            if key.lower() != "potcar":
    #                err = True
    #                logger.error(f"{key} not found")
    #            else:
    #                logger.warning(f"{key} not found")
    #    if err:
    #        sys.exit()

    #    flist = "\n".join([str(ii) for ii in results.values()])
    #    logger.debug(f"the following files were found:\n{flist}")

    return results


def load_poscar(value):
    if is_single_line(value):
        poscar = Poscar.from_file(
            value, check_for_potcar=False
        )  # set symbols manually later
    else:
        poscar = Poscar.from_str(
            value,
        )  # set symbols manually later
    return poscar


def load_incar(value):
    if is_single_line(value):
        incar = Incar.from_file(value)
    else:
        incar = Incar.from_dict(value)

    return incar


def load_kpoints(value):
    if is_single_line(value):
        kpoints = Kpoints.from_file(value)
    else:
        kpoints = Kpoints.from_str(value)
    return kpoints


def load_files(available_files):

    incar, poscar, potcar, kpoints = None, None, None, None
    if available_files["INCAR"] is not None:
        incar = load_incar(available_files["INCAR"])
    if available_files["POSCAR"] is not None:
        poscar = load_poscar(available_files["POSCAR"])

    if available_files["POTCAR"] is not None:
        potcar = Potcar.from_file(available_files["POTCAR"])
    if available_files["KPOINTS"] is not None:
        kpoints = load_kpoints(available_files["KPOINTS"])

    return incar, poscar, potcar, kpoints


def reorder_list(reference, target):
    mapping = {ref: i for i, ref in enumerate(reference)}
    return sorted(target, key=lambda x: mapping[x.split("_")[0]])


def prepare_potcar(cfg, symbols_poscar, potcar=None):
    potcar_sym = symbols_poscar
    functional = None

    if potcar is not None:
        functional = potcar.functional
        potcar_sym = potcar.symbols

    if cfg.calc.pseudo.variant is not None:
        potcar_sym = [str(ii) for ii in cfg.calc.pseudo.variant]
        elements_cfg = []
        for ii in potcar_sym:
            elements_cfg.append(ii.split("_")[0])

        if set(elements_cfg) != set(symbols_poscar):
            logger.error(
                f"poscar symbols {symbols_poscar} and pseudo variants {potcar_sym} are not compatible"
            )
            sys.exit()
        # for ss in symbols_poscar:
        #      potstr=''
        #      if len(str(cfg.calc.pseudo.variant)) >0:
        #            potstr=f"_{cfg.calc.pseudo.variant}"
        #      potcar_sym.append(f'{ss}{potstr}')

    if cfg.calc.functional is not None:
        functional = cfg.calc.functional
    potcar_sym = reorder_list(symbols_poscar, potcar_sym)

    potcar = Potcar(potcar_sym, functional=functional)
    sym_str = ""
    for ss in potcar_sym:
        sym_str += f" {ss}"

    logger.info(f"POTCAR generated: functional {functional}, symbols{sym_str}")

    return potcar


def compile_run_script(env, mpiexec, nproc, command, calcdir):
    lines = ["#!/bin/bash -l\n"]

    lines.append("\n" + env + "\n")
    cmd = f"{mpiexec} -np {nproc} {command}"

    lines.append(f"pushd {calcdir} || exit 1\n    " + cmd + "\npopd\n")

    return lines


def check_symbols_order(loop_result):
    for key, ii in loop_result.items():
        pot_sym = [ii.split("_")[0] for ii in ii.potcar.symbols]
        if ii.poscar.site_symbols != pot_sym:
            logger.error(
                f"symbols in POTCAR {pot_sym} and POSCAR {ii.poscar.site_symbols} are ordered differently"
            )
            sys.exit()


def compile_input_loop(
    cfg, incar, poscar, potcar, kpoints, file_poscar, parallelism_dict=None
):
    parameters_of_loops = []
    list_of_loops = []
    loop_result = {}
    name = cfg.dir.prefix
    messages = []

    if parallelism_dict is not None:
        incar.update(parallelism_dict)

    if cfg.loop is not None:
        name_tmp = name
        for ll in cfg.loop:
            if "list" in ll.interpolation:
                eoscar = Poscar.from_str("".join(poscar_lines))
                parameters_of_loops.append({"file": ll.file, "parameter": ll.parameter})
                values = ll.val
                list_of_loops.append(values)
            elif "interval" in ll.interpolation:
                parameters_of_loops.append({"file": ll.file, "parameter": ll.parameter})
                if "kpoints" in ll.parameter:
                    c_over_a = ll.get("c_over_a", 1)
                    logger.info(f"looping over KPOINTS with c/a={c_over_a:.2f}")
                    krange = np.arange(ll.val[0], ll.val[1], ll.val[2])
                    values = [[ii, ii, int(round(ii / c_over_a, 0))] for ii in krange]
                    if ll.get("include_gamma", True):
                        if not [1, 1, 1] in values:
                            values = [[1, 1, 1], *values]
                else:
                    values = np.arange(ll.val[0], ll.val[1], ll.val[2])
                list_of_loops.append(values)
            logger.info(f"looping on {ll.parameter.upper()}: {values}")

        list_of_calc = list(product(*list_of_loops))
        poscar_new = poscar
        for cc in list_of_calc:
            name_tmp = ""
            for idx, ii in enumerate(parameters_of_loops):
                dash = "" if idx == 0 else "-"
                if ii["file"].lower() == "incar":
                    name_tmp = (
                        name_tmp + f'{dash}{ii["parameter"].upper()}_{cc[idx]:.2f}'
                    )
                    incar.update({ii["parameter"].upper(): cc[idx]})
                if ii["file"].lower() == "poscar":
                    # special handling here for the lattice constant
                    name_tmp = name_tmp + f'{dash}{ii["parameter"]}_{cc[idx]:.2f}'
                    if ii["parameter"].lower() == "a":
                        with open(file_poscar, "r") as f:
                            poscar_lines = f.readlines()
                        poscar_lines[1] = f"{cc[idx]:.4f}\n"
                        poscar_new = Poscar.from_str("".join(poscar_lines))
                    if ii["parameter"].lower() == "volume":
                        poscar.structure.lattice = poscar.structure.lattice.scale(
                            cc[idx]
                        )
                        poscar_new = poscar
                    else:
                        err = NotImplementedError(
                            f'loop over parameter {ii["parameter"]} not implemented'
                        )
                        logger.error(err)
                        sys.exit()
                if ii["file"].lower() == "kpoints":
                    mesh_type = ii.get("mesh_type", "Gamma")
                    if "monk" in mesh_type.lower():
                        messages.append("KPOINTS mesh_type is Monkhorst-Pack")
                        kpoints = Kpoints.monkhorst_automatic(
                            tuple(cc[idx]),
                            comment="automatic Monkhorst-Pack Kpoint grid",
                        )
                    else:
                        messages.append("KPOINTS mesh_type is Gamma-centered")

                        kpoints = Kpoints.gamma_automatic(
                            tuple(cc[idx]),
                            comment="automatic Gamma-centered Kpoint grid",
                        )
                    name_tmp = (
                        name_tmp
                        + f"{dash}K{cc[idx][0]:d}_{cc[idx][1]:d}_{cc[idx][2]:d}"
                    )

            newinput = VaspInput(
                incar=copy.deepcopy(incar),
                poscar=copy.deepcopy(poscar_new),
                kpoints=copy.deepcopy(kpoints),
                potcar=copy.deepcopy(potcar),
            )
            loop_result.update({name_tmp: newinput})

    else:
        loop_result.update(
            {
                name: VaspInput(
                    incar=incar, poscar=poscar, kpoints=kpoints, potcar=potcar
                )
            }
        )

    check_symbols_order(loop_result)  # exit if symbols are not correctly set
    if len(messages) > 0:
        logger.info(messages[0])

    return loop_result


def join_if_list(inp, glue="-"):
    if inp is None:
        return inp
    elif isinstance(inp, str):
        return inp
    elif isinstance(inp, ListConfig):
        return glue.join([str(ii) for ii in inp if len(ii) > 0])


def set_directories(cfg, loop_result):
    prefix = join_if_list(cfg.dir.prefix)
    suffix = join_if_list(cfg.dir.suffix)
    rootdir = Path(os.getcwd()) / Path(prefix)
    subdir = "" if cfg.dir.subdir is None else join_if_list(cfg.dir.subdir)
    destname = suffix if suffix is not None else "./"
    workdir = rootdir / subdir / Path(destname)

    destinations = {}
    for name in loop_result.keys():
        destinations.update({name: workdir / Path(name)})
    return destinations


def write_calc(cfg, loop_result, destinations):
    for name, vaspinput in loop_result.items():
        calcdir = destinations[name]
        if decide_overwrite(calcdir, cfg.dir.overwrite):
            calcdir.mkdir(parents=True, exist_ok=True)
            logger.info(f"{calcdir}")
        else:
            logger.info("stopping")
            return 0

        vaspinput.write_input(calcdir)
    logger.info(f"parent directory: {calcdir.parent}")


def compile_sbatch_script(setup, env, srun_flags, command, name):

    if setup["job-name"] is None:
        setup["job-name"] = name[-12:]
    setup["error"] = str(name) + "/err"
    setup["output"] = str(name) + "/log"
    if setup["partition"] == "debug":
        setup["time"] = "00:30:00"

    sbatch_lines = ["#!/bin/bash -l\n"]

    for key, val in setup.items():
        if val is None:
            logger.warning(f"slurm setup {key} = {val}!")

        # deal with duplicate values (--hint can repeat)
        if not isinstance(val, str) and hasattr(val, "__len__") and (len(val) > 1):
            line = ""
            for ii in val:
                line += f"#SBATCH --{key}={ii}\n"
        else:
            line = f"#SBATCH --{key}={val}\n"
        sbatch_lines.append(line)

    sbatch_lines.append("\n" + env + "\n")

    # set srun flags

    result = ["srun"]
    for key, val in srun_flags.items():
        # deal with duplicate values
        if not isinstance(val, str) and hasattr(val, "__len__") and (len(val) > 1):
            for ii in val:
                result.append(f"--{key}={ii}")
        else:
            result.append(f"--{key}={val}")

    command = " ".join(result + [str(command)])
    sbatch_lines.append(f"pushd {name} || exit 1\n    " + command + "\npopd\n")
    sbatch_lines.append(f"\necho $SLURM_JOB_ID > run.{name}.success\n")

    return sbatch_lines


def write_exec_scripts(executor, loop_result, destinations, sbatch=False):
    for name, vaspinput in loop_result.items():
        calcdir = destinations[name]
        lines = None
        executable = False
        if OmegaConf.select(executor, "slurm") is not None:
            settings = executor.slurm
            lines = compile_sbatch_script(
                copy.deepcopy(settings.setup),
                settings.env,
                settings.srun_flags,
                settings.cmd,
                calcdir.name,
            )
        elif OmegaConf.select(executor, "local") is not None:
            settings = executor.local
            lines = compile_run_script(
                env=settings.env,
                mpiexec=settings.mpiexec,
                nproc=settings.nproc,
                command=settings.cmd,
                calcdir=calcdir.name,
            )
            executable = True
        file = calcdir.parent / Path(f"run.{name}")
        with open(file, "w") as f:
            f.writelines(lines)
        if executable:
            os.chmod(
                file, os.stat(file).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
            )
        if sbatch:
            script_dir = calcdir.parent
            result = subprocess.run(["sbatch", os.path.basename(file)], cwd=script_dir, capture_output=True, text=True)
            if result.returncode == 0:
                logger.info(f"Job submitted successfully: {result.stdout.strip()}")
            else:
                logger.error(f"Error submitting job: {result.stderr.strip()}")



def write_dpdisp_tasks(executor, loop_result, destinations):
    tasks = []
    for name, vaspinput in loop_result.items():
        calcdir = destinations[name]
        tt = {
            "command": "mpirun -np ${DPDISPATCHER_CPU_PER_NODE} vasp_std",
            "task_work_path": str(calcdir),
        }
        tasks.append(tt)
    file = calcdir.parent / Path(f"tasks.yaml")
    with open(file, "w") as f:
        yaml.dump(tasks, f, default_flow_style=False, indent=4)
        # json.dump(tasks, f, indent=4)


def set_parallelism_dict(cfg):
    para_dict = None
    if cfg.executor is not None:
        if "daint" in cfg.executor:
            exe = cfg.executors[cfg.executor]
            ngpus = exe.slurm.setup.nodes * 4
            para_dict = {"KPAR": ngpus}

    return para_dict


@hydra.main(config_path="config", config_name="config", version_base=None)
def main(cfg):

    logger.info("startup")

    # scan for available files
    logger.debug(
        f'scanning for input {" ".join([ ii for ii in cfg.dir if ii is not None])}'
    )
    available_files = scan_vasp_files(
        ["INCAR", "POSCAR", "KPOINTS", "POTCAR"], cfg.source
    )

    # load parts of the input
    incar, poscar, potcar, kpoints = load_files(available_files)
    incar_cfg = cfg.get("incar", None)
    if incar_cfg is not None:
        incar = Incar(incar_cfg)
        msg = "; ".join(str(incar).split("\n"))
        logger.info("loading INCAR from config file:")
        logger.info(f"{msg}")
    else:
        msg = "; ".join(str(incar).split("\n"))
        logger.info(f"loading INCAR:")
        logger.info(f"{msg}")

    # prepare potcar, from config spec or from input
    potcar = prepare_potcar(cfg, poscar.site_symbols, potcar)

    # prepare parallelization flags for incar file
    para_dict = set_parallelism_dict(cfg)

    # loop over parameters
    loop_result = compile_input_loop(
        cfg,
        incar,
        poscar,
        potcar,
        kpoints,
        available_files["POSCAR"],
        parallelism_dict=para_dict,
    )

    # check poscar/potcar species order

    # set directories
    destinations = set_directories(cfg, loop_result)

    # write inputs
    write_calc(cfg, loop_result, destinations)

    # write slurm scripts
    if cfg.executor is not None:
        exe = cfg.executors[cfg.executor]
        sbatch=cfg.get('sbatch',False)
        write_exec_scripts(exe, loop_result, destinations, sbatch)
        if "csea" in cfg.executor:
            write_dpdisp_tasks(exe, loop_result, destinations)
    else:
        logger.warning("no executor specified - unable to write run scripts")

    plural = "" if len(loop_result) < 2 else "s"
    logger.info(f"{len(loop_result)} folder{plural} total")


# da fare ora:

# implementare help message
# clean exif from 'do you want to continue'
# test config: subdir, prefix, source.*, executor
# TODO: fix loop over lattice constant when poscar is from config file
# add print config.yaml to target directory
# add feature recreate from old config file
# make incar keys always uppercase


# DONE implementare precisione variabile per i parametri di loop
# DONE imparare come impostare il loop dalla cli
# DONE sistemare config file
# DONE scrivere executor per slurm
# DONE scrivere executor per cseasrv
# DONE sistemare root dir
# DONE creare il naming system e la cartella
# DONE caricare files su pymatgen
# DONE scrivere files
# DONE creare un sistema per caricare files da esterno
# DONE (potrebbe essere un sistema che prenda un template e lo scriva sul posto, tipo py vasp.py write-template=True, scrive la cartella, poi io modifico l'incar, poi runno py vasp.py source.dir=cartella_appena_modificata source.poscar=eccetera)


if __name__ == "__main__":
    main()
