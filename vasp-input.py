import hydra
from pathlib import Path
from pymatgen.io.vasp import Poscar, Kpoints, Incar, Potcar, VaspInput
from omegaconf import OmegaConf
from omegaconf.listconfig import ListConfig
import numpy as np
import sys
from itertools import product
import json
import copy
import os, stat

from logger_setup import logger




def ask_if_overwrite(dd):
      folder=Path(dd)
      if folder.exists():
            user_input = input(f"{str(folder)} exists. Do you want to continue? (y/n): ")
            if user_input.lower() == 'y':
              proceed=True
            else:
              proceed=False
      else:
            proceed = True
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
            logger.warning('source.dir not set')
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
                  logger.error(f'{key} not found')
                  if key != 'potcar':
                        err=True
      if err:
            sys.exit()
      
      flist="\n".join([ str(ii) for ii in results.values()])
      logger.debug(f'the following files were found:\n{flist}')
      
      return results

def load_files(available_files):
      
      incar=Incar.from_file(available_files['INCAR'])
      poscar=Poscar.from_file(available_files['POSCAR'])
      potcar=Potcar.from_file(available_files['POTCAR'])
      kpoints=Kpoints.from_file(available_files['KPOINTS'])

      return incar, poscar, potcar, kpoints

def prepare_potcar(cfg, symbols, potcar=None):
      potcar_sym=[]
      functional=None
      if potcar is not None:
            functional=potcar.functional
            potcar_sym=potcar.symbols

      if cfg.calc.pseudo.variant is not None:
            potcar_sym=[]
            for ss in symbols:
                  potstr=''
                  if len(str(cfg.calc.pseudo.variant)) >0:
                        potstr=f"_{cfg.calc.pseudo.variant}"
                  potcar_sym.append(f'{ss}{potstr}')

      if cfg.calc.functional is not None:
            functional=cfg.calc.functional

      potcar=Potcar(potcar_sym, functional=functional)
      sym_str=''
      for ss in potcar_sym:
            sym_str += f' {ss}'
      
      logger.info(f'POTCAR generated: functional {functional}, symbols{sym_str}')

      return potcar





def compile_run_script(env, mpiexec, nproc, command, calcdir):
      lines=['#!/bin/bash -l\n']

      
      lines.append('\n'+env+'\n')
      cmd = f'{mpiexec} -np {nproc} {command}'
      
      lines.append(f'pushd {calcdir} || exit 1\n    ' + cmd + '\npopd\n')
      
      return lines


def compile_input_loop(cfg, incar, poscar, potcar, kpoints, file_poscar):
      parameters_of_loops=[]
      list_of_loops=[]
      loop_result={}
      name = cfg.dir.prefix
      
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
                  logger.info(f'looping on {ll.parameter.upper()}: {values}')
      
            list_of_calc=list(product(*list_of_loops))
            for cc in list_of_calc:
                  name_tmp=''      
                  for idx, ii in enumerate(parameters_of_loops):
                        dash='' if idx == 0 else '-'
                        if ii['file'].lower() == 'incar':
                              name_tmp=name_tmp+f'{dash}{ii["parameter"].upper()}_{cc[idx]:.2f}'
                              incar.update({ii["parameter"].upper():cc[idx]})
                        if ii['file'].lower() == 'poscar':
                              # special handling here for the lattice constant
                              name_tmp=name_tmp+f'{dash}{ii["parameter"]}_{cc[idx]:.2f}' 
                              if ii["parameter"].lower() == 'a':
                                    with open(file_poscar,'r') as f:
                                          poscar_lines=f.readlines()
                                    poscar_lines[1]=f'{cc[idx]:.4f}\n'
                              else:
                                    err=NotImplementedError(f'loop over parameter {ii["parameter"]} not implemented')
                                    logger.error(err)
                                    sys.exit()
                        poscar=Poscar.from_str("".join(poscar_lines))
                              

                  newinput=VaspInput(incar=copy.deepcopy(incar),poscar=copy.deepcopy(poscar), \
                                                kpoints=copy.deepcopy(kpoints), potcar=copy.deepcopy(potcar))                 
                  loop_result.update({name_tmp:newinput})

      else:
            loop_result.update({name:VaspInput(incar=incar,poscar=poscar, kpoints=kpoints, potcar=potcar)})
      
      return loop_result

def join_if_list(inp, glue='-'):
      if inp is None:
            return inp
      elif isinstance(inp, str):
            return inp
      elif isinstance(inp, ListConfig):
             return glue.join([str(ii) for ii in inp if len(ii)>0])

def set_directories(cfg, loop_result):
      prefix=join_if_list(cfg.dir.prefix)
      suffix=join_if_list(cfg.dir.suffix)
      rootdir = Path(os.getcwd())/Path(prefix) 
      subdir = '' if cfg.dir.subdir is None else join_if_list(cfg.dir.subdir)
      destname=suffix if suffix is not None else './'
      workdir= rootdir/subdir/Path(destname)
      
      destinations={}
      for name in loop_result.keys():
            destinations.update({name:workdir/Path(name)})
      return destinations

def write_calc(cfg, loop_result, destinations):
      for name, vaspinput in loop_result.items():
            calcdir=destinations[name]      
            if decide_overwrite(calcdir, cfg.dir.overwrite):
                  calcdir.mkdir(parents=True, exist_ok=True)
                  logger.info(f'{calcdir}')
            else:
                  logger.info('stopping')
                  return 0
            
            vaspinput.write_input(calcdir)

def compile_sbatch_script(setup, env, command,  name):

      if setup['job-name'] is None:
            setup['job-name']=name[-12:]
      setup['error']=str(name)+'/err'
      setup['output']=str(name)+'/log'
      if setup['partition'] == 'debug':
            setup['time']='00:30:00'
        
      sbatch_lines=['#!/bin/bash -l\n']

      for key, val in setup.items():
            if val is None:
                  logger.warning(f'slurm setup {key} = {val}!')

            # deal with duplicate values (--hint can repeat)
            if not isinstance(val, str) and hasattr(val, "__len__") and (len(val) > 1):
                  line=''
                  for ii in val:
                        line+=f'#SBATCH --{key}={ii}\n'
            else:
                  line=f'#SBATCH --{key}={val}\n'
            sbatch_lines.append(line)
      
      sbatch_lines.append('\n'+env+'\n')
      
      sbatch_lines.append(f'pushd {name} || exit 1\n    ' + command + '\npopd\n')
      sbatch_lines.append(f'\necho $SLURM_JOB_ID > run.{name}.success\n')
      
      return sbatch_lines

def write_exec_scripts(executor, loop_result, destinations):
    for name, vaspinput in loop_result.items():
          calcdir=destinations[name]
          lines = None
          executable=False
          if OmegaConf.select(executor, "slurm") is not None:
              settings=executor.slurm
              lines=compile_sbatch_script(copy.deepcopy(settings.setup), settings.env, settings.cmd, calcdir.name)
          elif OmegaConf.select(executor, "local") is not None:
              settings=executor.local
              lines=compile_run_script(env=settings.env, mpiexec=settings.mpiexec, nproc=settings.nproc, command=settings.cmd, \
                      calcdir=calcdir.name)
              executable=True
          file=calcdir.parent/Path(f'run.{name}')
          with open(file, 'w') as f:
                f.writelines(lines)
          if executable:
              os.chmod(file, os.stat(file).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)



@hydra.main(config_path="config", config_name="config", version_base=None)
def main(cfg):

      logger.info('startup')
     
      # scan for available files 
      logger.debug(f'scanning for input {" ".join([ ii for ii in cfg.source.values() if ii is not None])}')
      available_files=scan_vasp_files(['INCAR', 'POSCAR', 'KPOINTS','POTCAR'], cfg.source)

      # load parts of the input
      incar, poscar, potcar, kpoints = load_files(available_files)
      logger.debug(f'the reference INCAR is:\n{str(incar)[:-1]}')

      # prepare potcar, from config spec or from input
      potcar = prepare_potcar(cfg, poscar.structure.symbol_set, potcar) 
      
      # loop over parameters
      loop_result = compile_input_loop(cfg, incar, poscar, potcar, kpoints, available_files['POSCAR'])

      # set directories
      destinations=set_directories(cfg, loop_result)

      # write inputs
      write_calc(cfg, loop_result, destinations)            
      
      # write slurm scripts
      if cfg.executor is not None:
            exe=cfg.executors[cfg.executor]
            write_exec_scripts(exe, loop_result, destinations)
      else:
            logger.warning('no executor specified - unable to write run scripts')


      plural='' if len(loop_result) < 2 else 's'
      logger.info(f'{len(loop_result)} folder{plural} total')
     
      




            
      

# da fare ora:

# implementare help message
# clean exif from 'do you want to continue'
# test config: subdir, prefix, source.*, executor


 
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
      


if __name__ == '__main__':
      main()
