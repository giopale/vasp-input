import hydra
from pathlib import Path
from pymatgen.io.vasp import Poscar, Kpoints, Incar, Potcar, VaspInput
import numpy as np
import sys
from itertools import product
import json
import copy


import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,  # Set the logging level
    format="%(asctime)s - %(levelname)s - %(message)s",  # Log format
    datefmt="%Y-%m-%d %H:%M:%S",  # Date format
    handlers=[
        logging.FileHandler("vasp-input.log"),  # Log to a file
        logging.StreamHandler()  # Log to console
    ]
)



def ask_if_overwrite(dd):
      folder=Path(dd)
      user_input = input(f"{folder.name} exists. Do you want to continue? (y/n): ")
      if user_input.lower() == 'y':
        proceed=True
      else:
        proceed=False
      return proceed

def decide_overwrite(dd,overwrite):
      if overwrite:
            return True
      else:
            return ask_if_overwrite(dd)


def scan_vasp_files(keys,dict_files):
      # scan a directory for desired files
      # if specified, find other files and add them (with priority)


      if dict_files.dir is None:
            logging.warning('source.dir not set')
      sourcedir = Path(dict_files.dir)
      results={}

      for key in keys:
            value=next(sourcedir.glob(f"{key}*"), None)
            results.update({key:value})
      
      for key, value in dict_files.items(): # <--------------------------- to test
            if key == 'dir':
                  continue
            elif value is not None:
                  if Path(value).exists():
                        available_files.update({key.upper():Path(value)})
      err=False
      for key,val in results.items():
            if val is None:
                  logging.error(f'{key} not found')
                  if key != 'potcar':
                        err=True
      if err:
            sys.exit()
      
      flist="\n".join([ str(ii) for ii in results.values()])
      logging.debug(f'the following files were found:\n{flist}')
      
      return results

def load_files(available_files):
      
      incar=Incar.from_file(available_files['INCAR'])
      poscar=Poscar.from_file(available_files['POSCAR'])
      potcar=Potcar.from_file(available_files['POTCAR'])
      kpoints=Kpoints.from_file(available_files['KPOINTS'])

      return incar, poscar, potcar, kpoints

def prepare_potcar(cfg, symbols, potcar=None):
      potcar_sym=[]
      for ss in symbols:
            potstr=''
            if len(str(cfg.pseudo.variant)) >0:
                  potstr=f"_{cfg.pseudo.variant}"
            potcar_sym.append(f'{ss}{potstr}')

      if potcar is None or cfg.pseudo.variant is not None:
            potcar=Potcar(potcar_sym, functional=cfg.functional)
            sym_str=''
            for ss in potcar_sym:
                  sym_str += f' {ss}'
            
            logging.info(f'POTCAR {cfg.functional} generated:{sym_str}')
      else:
            sym_str=''
            for ss in potcar.symbols:
                  sym_str += f' {ss}'
            logging.info(f'POTCAR {potcar.functional} from input{sym_str}')


      return potcar



def compile_sbatch_script(setup, env, command, calcdir):
      sbatch_lines=['#!/bin/bash -l\n']

      for key, val in setup.items():
            if val is None:
                  logging.warning(f'slurm setup {key} = {val}!')

            # deal with duplicate values (--hint can repeat)
            if not isinstance(val, str) and hasattr(val, "__len__") and (len(val) > 1):
                  line=''
                  for ii in val:
                        line+=f'#SBATCH --{key}={ii}\n'
            else:
                  line=f'#SBATCH --{key}={val}\n'
            sbatch_lines.append(line)
      
      sbatch_lines.append('\n'+env+'\n')
      
      sbatch_lines.append(f'pushd {calcdir} || exit 1\n    ' + command + '\npopd\n')
      
      return sbatch_lines

def compile_input_loop(cfg, incar, poscar, potcar, kpoints, file_poscar):
      parameters_of_loops=[]
      list_of_loops=[]
      loop_result={}
      name = cfg.prefix
      
      if cfg.loop is not None:
            name_tmp=name
            for ll in cfg.loop:
                  if 'list' in ll.interpolation:
                        parameters_of_loops.append({'file':ll.file,'parameter':ll.parameter})
                        values=ll.val
                        list_of_loops.append(values)
                  elif 'interval' in ll.interpolation:
                        parameters_of_loops.append({'file':ll.file,'parameter':ll.parameter})
                        values=np.arange(ll.val[0], ll.val[1], ll.val[2])
                        list_of_loops.append(values)
                  logging.info(f'looping on {ll.parameter.upper()}: {values}')
      
            list_of_calc=list(product(*list_of_loops))
            for cc in list_of_calc:
                  name_tmp=''      
                  for idx, ii in enumerate(parameters_of_loops):
                        dash='' if idx == 0 else '-'
                        if ii['file'].lower() == 'incar':
                              name_tmp=name_tmp+f'{dash}{ii["parameter"].upper()}_{cc[idx]:.2g}'
                              incar.update({ii["parameter"].upper():cc[idx]})
                        if ii['file'].lower() == 'poscar':
                              # special handling here for the lattice constant
                              name_tmp=name_tmp+f'{dash}{ii["parameter"]}_{cc[idx]:.2g}'
                              if ii["parameter"].lower() == 'a':
                                    with open(file_poscar,'r') as f:
                                          poscar_lines=f.readlines()
                                    poscar_lines[1]=f'{cc[idx]:.4f}\n'
                              else:
                                    err=NotImplementedError(f'loop over parameter {ii["parameter"]} not implemented')
                                    logging.error(err)
                                    sys.exit()
                        poscar=Poscar.from_str("".join(poscar_lines))
                              

                  newinput=VaspInput(incar=copy.deepcopy(incar),poscar=copy.deepcopy(poscar), \
                                                kpoints=copy.deepcopy(kpoints), potcar=copy.deepcopy(potcar))                 
                  loop_result.update({name_tmp:newinput})

      else:
            loop_result.update({name:VaspInput(incar=incar,poscar=poscar, kpoints=kpoints, potcar=potcar)})
      
      return loop_result

def set_directories(cfg, loop_result):
      suffix='' if cfg.suffix is None else '-'+cfg.suffix
      rootdir = Path(cfg.rootdir)/Path(cfg.prefix+suffix) 
      subdir = '' if cfg.subdir is None else cfg.subdir
      workdir= rootdir/subdir
      
      destinations={}
      for name in loop_result.keys():
            destinations.update({name:workdir/Path(name)})
      return destinations

def write_calc(cfg, loop_result, destinations):
      for name, vaspinput in loop_result.items():
            calcdir=destinations[name]      
            if decide_overwrite(calcdir, cfg.overwrite):
                  calcdir.mkdir(parents=True, exist_ok=True)
                  logging.info(f'{calcdir}')
            else:
                  logging.info('stopping')
                  return 0
            
            vaspinput.write_input(calcdir)

def write_slurm_scripts(cfg, loop_result, destinations):
      for name, vaspinput in loop_result.items():
            calcdir=destinations[name]
            sbatchconfig=copy.deepcopy(cfg.slurm.setup)
            if sbatchconfig['job-name'] is None:
                  sbatchconfig['job-name']=name[-12:]
            sbatchconfig['error']=str(calcdir.name)+'/err'
            sbatchconfig['output']=str(calcdir.name)+'/log'
            if sbatchconfig['partition'] == 'debug':
                  sbatchconfig['time']='00:30:00'

            sbatch_lines=compile_sbatch_script(sbatchconfig, cfg.slurm.env, cfg.slurm.cmd, calcdir)
            sbatch_lines.append(f'\necho $SLURM_JOB_ID > run.{name}.success\n')
            sbatch_file=calcdir.parent/Path(f'run.{name}')
            with open(sbatch_file, 'w') as f:
                  f.writelines(sbatch_lines)



@hydra.main(config_path="config", config_name="config", version_base=None)
def main(cfg):

      logging.info('startup')
     
      # scan for available files 
      logging.debug(f'scanning for input {" ".join([ ii for ii in cfg.source.values() if ii is not None])}')
      available_files=scan_vasp_files(['INCAR', 'POSCAR', 'KPOINTS','POTCAR'], cfg.source)

      # load parts of the input
      incar, poscar, potcar, kpoints = load_files(available_files)
      logging.debug(f'the reference INCAR is:\n{str(incar)[:-1]}')

      # prepare potcar, from config spec or from input
      potcar = prepare_potcar(cfg, poscar.structure.symbol_set, potcar) 
      
      # loop over parameters
      loop_result = compile_input_loop(cfg, incar, poscar, potcar, kpoints, available_files['POSCAR'])

      # set directories
      destinations=set_directories(cfg, loop_result)

      # write inputs
      write_calc(cfg, loop_result, destinations)            
      
      # write slurm scripts
      write_slurm_scripts(cfg, loop_result, destinations)

      plural='' if len(loop_result) < 2 else 's'
      logging.info(f'{len(loop_result)} folder{plural} total')
     
      




            
      

 # da fare ora:
# DONE creare il naming system e la cartella
# DONE caricare files su pymatgen
# DONE scrivere files
# DONE creare un sistema per caricare files da esterno
# DONE (potrebbe essere un sistema che prenda un template e lo scriva sul posto, tipo py vasp.py write-template=True, scrive la cartella, poi io modifico l'incar, poi runno py vasp.py source.dir=cartella_appena_modificata source.poscar=eccetera)
      


if __name__ == '__main__':
      main()
