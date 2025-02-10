# vasp-input
A python-based package to generate VASP input files.


## Features

- parse template input files `INCAR`, `KPOINTS`, `POSCAR`, and `POTCAR`
- loop over `INCAR` parameters
- loop over lattice parameter in `POSCAR`
- generate `slurm` scripts or `bash` scripts to run the calculations

## Installation

This package requires `hydra` for building the configuration and `pymatgen` to manage the vasp input files.

To install `hydra` use ```pip install hydra-core```.

Installing pymatgen requires two steps:
1) Install the `pymatgen` package with `pip install pymatgen`
2) Set up the pseudopotentials: after the installation run ``` pmg config -p <EXTRACTED_VASP_POTCAR> <MY_PSP>```where `<EXTRACTED_VASP_POTCAR>` is the path to the extracted VASP pseudopotential files as obtained from VASP, and `<MY_PSP>` is the desired path where you would like to store the reformatted, `pymatgen`-compatible pseudopotential files. See the [pymatgen website](https://pymatgen.org/installation.html#potcar-setup) for more.


## Usage

### Command line interface
Set up a template folder with the vasp input files:
```
template/
├── INCAR
├── KPOINTS
├── POSCAR
└── POTCAR
```

then, from the command line, run:

```bash
python3 <path_to_source>/vasp-input.py \
    source.dir=template \
    dir.prefix=Si \
    dir.subdir=PBE \
    dir.suffix=cell-search \
    loop='[{file:poscar,parameter:a,interpolation:interval,val:[5.5,5.75,0.1]}]'  \ #python syntax for dictionary
    dir.overwrite=True \
    executor=daint-gh200
```
You'll obtain the following tree:

```
Si
└── PBE
    └── cell-search
        ├── a_5.50
        │   ├── INCAR
        │   ├── KPOINTS
        │   ├── POSCAR
        │   └── POTCAR
        ├── a_5.60
        │   ├── INCAR
        │   ├── KPOINTS
        │   ├── POSCAR
        │   └── POTCAR
        ├── a_5.70
        │   ├── INCAR
        │   ├── KPOINTS
        │   ├── POSCAR
        │   └── POTCAR
        ├── run.a_5.50
        ├── run.a_5.60
        └── run.a_5.70
```


### Configuration file

Alternatively, you can create a configuration file `config.yaml`:

```yaml

source:
  dir: template
  poscar: null
  kpoints: null
  incar: null
  potcar: null
dir:
  prefix: Si
  suffix: cell-search
  subdir: PBE
  overwrite: true
loop:
- file: poscar
  parameter: a
  interpolation: interval
  val:
  - 5.5
  - 5.75
  - 0.1
calc:
  functional: PBE
  pseudo:
    variant: null
executor: daint-gh200

executors:
  daint-gh200:
    slurm:
      cmd: srun --cpu-bind=socket ~/mps-wrapper.sh vasp_std
      setup:
        job-name: null
        account: lp07
        constraint: gpu
        hint:
        - nomultithread
        - exclusive
        nodes: 1
        ntasks-per-node: 32
        ntasks-per-core: 1
        cpus-per-task: 1
        partition: normal
        time: '12:00:00'
        output: log
        error: err
        uenv: vasp/v6.4.3:v1
      env: |
        export PATH=/user-environment/env/vasp/bin/:$PATH
        export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
        export MPICH_MALLOC_FALLBACK=1
        export MPICH_GPU_SUPPORT_ENABLED=1

```

then run:

```bash
python3 <path_to_source>/vasp-input.py --config-path=$(pwd) --cofig-name=config.yaml
```




## Configuration

The program is configured using a `yaml` file. You can print a template of the basic config file with:

```bash
python3 vasp-input.py --cfg=job
```

The output is going to look like:

```yaml
executors:
  daint-gh200:
    # ... see later
  cseasrv:
    # ... see later
  
source: # path to the template files
  dir: null
  poscar: null
  kpoints: null
  incar: null
  potcar: null
dir: # output directory settings
  prefix: null
  suffix: null
  subdir: null
  overwrite: false
loop: # loop settings
    - file: null
      parameter: null
      interpolation: null
      val: null
calc: # calculation settings
  functional: PBE
  pseudo:
    variant: null
```


The configuration file is in `yaml` format.
The settings are divided into sections:
- `executors` contains the settings for the different clusters where the calculations are going to be run.
- `source` contains the path to the template files.
- `dir` contains the settings for the output directories.
- `loop` contains the settings for the loops over the lattice parameter in `POSCAR` and/or the parameters in `INCAR`.
- `calc` contains the settings for the calculation.

### Template files
The program runs on a set of template files. You can provide either the path of a directory containing the template files or the path of a single file.
Single files have priority over the full directory.

for example, the following configuration:

```yaml
source:
  dir: /home/<user>/VASP/Si-origin
  poscar: /home/<user>/VASP/structures/Si/POSCAR
  kpoints: null
  incar: null
  potcar: null
```
will cause the program to scan `INCAR`, `KPOINTS`, `POSCAR`, and `POTCAR` files in `/home/<user>/VASP/Si-origin` and use the `POSCAR` file in `/home/<user>/VASP/structures/Si/POSCAR`.


### Output directories
The program will create a directory for each calculation.
The directory name is composed of the `prefix`, the `suffix`, and the `subdir` settings.
The `prefix` and `suffix` are strings that are added to the directory name.
The `subdir` is a string that is added to the directory name and is used to create a subdirectory.
The `overwrite` setting is a boolean that specifies if the program should overwrite the directory if it already exists.

For example, the following configuration:

```yaml
dir:
  prefix: Si
  subdir: PBE
  suffix: cell-search
  overwrite: false
```

will produce the following tree:
```
Si
└── PBE
    └── cell_search
        ├── ...
```

Only **prefix is mandatory.**
If no `subdir` or `suffix` are provided, the calculation directories will be created in the `prefix` directory.

### Calculation
The program re-generates a POTCAR file for each calculation, using `pymatgen`.
```yaml
calc:
  functional: PBE
  pseudo:
    variant: null
```
With `variant` can choose between the different variants of the pseudopotentials provided by vasp, e.g. `sv_GW`.
Set to `null` to use the variant provided in the template POTCAR file parsed at the beginning of the program.

### Loop over parameters

The program can loop over the lattice parameter in the `POSCAR` file and/or the parameters in the `INCAR` file.

At the moment, two kinds of loops are supported:

- loop over lattice parameter 'a' in `POSCAR` file
- loop over any parameter in `INCAR` file

To set up a loop for the lattice constant:
```yaml
loop:
    - file: poscar
      parameter: a
      interpolation: interval
      val: [5.5, 5.7, 0.05]
```

To set up a loop for the parameters in `INCAR`:

```yaml
loop:
    - file: incar
      parameter: ENCUT
      interpolation: list
      val: [50, 100, 200, 400]
```


Two kinds of `interpolation` are supported:
- `interval` for a linear interpolation between the values in `val`; the first value is the starting point, the second value is the end point, and the third value is the step (see `numpy.arange`)
- `list` for a list of values in `val`

Multiple loops can be set combined in the `loop` section. The program will loop over all the combinations of the parameters and name the folders accordingly.

```yaml
loop:
    - file: incar
      parameter: ENCUT
      interpolation: list
      val: [50, 100, 200, 400]
    - file: poscar
      parameter: a
      interpolation: interval
      val: [5.5, 5.7, 0.05]
```


### Executors

The program can generate `slurm` scripts or `bash` scripts to run the calculations.
The `executors` section contains the settings for the different clusters where the calculations are going to be run.

There are two types of executors:

- `slurm` for the `slurm` scheduler. The program will generate a `slurm` script for each calculation.
- `local` for the local machine. The program will generate a `bash` script for each calculation.

For the slurm executor, the current template is:
```yaml
executors:
      daint-gh200:
            slurm:
                  cmd: srun --cpu-bind=socket ~/mps-wrapper.sh vasp_std
                  setup:
                        job-name: null
                        account: lp07
                        constraint: gpu
                        hint: [nomultithread, exclusive]
                        nodes: 1
                        ntasks-per-node: 32
                        ntasks-per-core: 1
                        cpus-per-task: 1
                        partition: normal
                        time: '12:00:00'
                        output: log
                        error: err
                        uenv: 'vasp/v6.4.3:v1'

                  env: |
                      export PATH=/user-environment/env/vasp/bin/:$PATH
                      export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
                      export MPICH_MALLOC_FALLBACK=1
                      export MPICH_GPU_SUPPORT_ENABLED=1
```

Further parameters can be added to the `slurm.setup` section.

The resulting slurm script is:

```bash
#!/bin/bash -l
#SBATCH --job-name=a_5.30
#SBATCH --account=lp07
#SBATCH --constraint=gpu
#SBATCH --hint=nomultithread
#SBATCH --hint=exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --time=12:00:00
#SBATCH --output=a_5.30/log
#SBATCH --error=a_5.30/err
#SBATCH --uenv=vasp/v6.4.3:v1

export PATH=/user-environment/env/vasp/bin/:$PATH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MPICH_MALLOC_FALLBACK=1
export MPICH_GPU_SUPPORT_ENABLED=1

pushd /capstor/scratch/cscs/pgiorgio/VASP/Si/cell-search/a_5.30 || exit 1
    srun --cpu-bind=socket ~/mps-wrapper.sh vasp_std
popd

echo $SLURM_JOB_ID > run.a_5.30.success
```

For running on the local machine, the current template is:
```yaml
executors:
      cseasrv:
            local:
                  cmd: vasp_std
                  nproc: 128
                  mpiexec: mpirun
                  env: |
                      export SPACK_ROOT=/home/palermo/spack;
                      export SPACK_COLOR=always;
                      source \${SPACK_ROOT}/share/spack/setup-env.sh
                      spack load /cptwhgb #vasp@6.4.3~cuda+fftlib+hdf5+openmp+shmem
```
which corresponds to this script:

```bash
#!/bin/bash -l

export SPACK_ROOT=/home/palermo/spack;
export SPACK_COLOR=always;
source ${SPACK_ROOT}/share/spack/setup-env.sh
spack load /cptwhgb #vasp@6.4.3~cuda+fftlib+hdf5+openmp+shmem

pushd /home/palermo/VASP/Si/cell_search/a_5.41 || exit 1
    mpirun -np 128 vasp_std
popd
```

You can add a new executor to the `executors` dictionary. Use the templates provided for either the `slurm` or `local` execution.






















