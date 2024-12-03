![MRChem logo](https://github.com/MRChemSoft/mrchem/raw/master/doc/gfx/logo_full.png)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3606658.svg)](https://doi.org/10.5281/zenodo.3606658)
[![License](https://img.shields.io/badge/license-%20LGPLv3-blue.svg)](../master/LICENSE)
[![Documentation Status](https://readthedocs.org/projects/mrchem/badge/?version=latest)](http://mrchem.readthedocs.io/en/latest/?badge=latest)
![Build and test MRChem](https://github.com/MRChemSoft/mrchem/workflows/Build%20and%20test%20MRChem/badge.svg)
[![CircleCI](https://circleci.com/gh/MRChemSoft/mrchem/tree/master.svg?style=svg)](https://circleci.com/gh/MRChemSoft/mrchem/tree/master)
[![codecov](https://codecov.io/gh/MRChemSoft/mrchem/branch/master/graph/badge.svg)](https://codecov.io/gh/MRChemSoft/mrchem)

MRChem is a numerical real-space code for molecular electronic structure
calculations within the self-consistent field (SCF) approximations of quantum
chemistry (Hartree-Fock and Density Functional Theory).

The code is being developed at the Hylleraas Centre for Quantum Molecular
Sciences at UiT - The Arctic University of Norway.

### Documentation: [mrchem.readthedocs.io](http://mrchem.readthedocs.io)


## Installation

For optimal performance it is recommended to build from source, as the packaged
builds are quite generic without architecture specific optimizations.


### From source

To build MRChem from source with MPI+OpenMP parallelization:

    $ git clone https://github.com/MRChemSoft/mrchem.git
    $ cd mrchem
    $ ./setup --prefix=<install-dir> --omp --mpi --cxx=<mpi-compiler> <build-dir>
    $ cd <build-dir>
    $ make
    $ make test
    $ make install

All dependencies will be fetched at configure time, if not already available.
For more information on different kinds of builds, see
[installation instructions](http://mrchem.readthedocs.io/en/latest/installation.html).


### Using Conda

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/mrchem/badges/version.svg)](https://anaconda.org/conda-forge/mrchem)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/mrchem/badges/latest_release_date.svg)](https://anaconda.org/conda-forge/mrchem)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/mrchem/badges/downloads.svg)](https://anaconda.org/conda-forge/mrchem)

To install MRChem in a Conda environment `myenv`:

    $ conda create -n myenv
    $ conda activate myenv
    $ conda install -c conda-forge mrchem               # latest version (OpenMP only)
    $ conda install -c conda-forge mrchem=1.0.0         # tagged version (OpenMP only)
    $ conda install -c conda-forge mrchem=*=*openmpi*   # latest version (MPI+OpenMP)
    $ conda install -c conda-forge mrchem=*=*mpich*     # latest version (MPI+OpenMP)

To list all available versions

    $ conda search -c conda-forge mrchem


### Using Spack

To install MRChem in a Spack environment `myenv`:

    $ spack env create myenv
    $ spack env activate myenv
    $ spack install mrchem                              # latest version (MPI+OpenMP)
    $ spack install mrchem @1.0.0                       # tagged version (MPI+OpenMP)
    $ spack install mrchem -mpi                         # latest version (OpenMP only)

For information on available Spack builds:

    $ spack info mrchem


### Using EasyBuild

To install MRChem in an EasyBuild/Lmod environment (only MPI+OpenMP version
available):

    $ eb MRChem-<version>-<toolchain> --fetch
    $ eb MRChem-<version>-<toolchain> --robot
    $ module load MRChem/<version>-<toolchain>

See
[EasyBuild](https://github.com/easybuilders/easybuild-easyconfigs/tree/develop/easybuild/easyconfigs/m/MRChem)
for available `<versions>` and `<toolchains>`.


### Using Singularity

Singularity recipe files are provided under `recipes/` for building local container images using
the current state of the source. Requires Singularity >= v3.2 as well as `sudo` rights on the
machine you are building on:

    $ sudo singularity build <image_name>.sif recipes/Singularity.<variant>

Recipes are provided for a pure OpenMP build (`recipes/Singularity.nompi`) and one MPI+OpenMP version,
using `OpenMPI-4.0` (`recipes/Singularity.openmpi4.0`).

Official MRChem images can also be downloaded from the GitHub Container Registry.

Latest `master` version (here OpenMP variant):

    $ singularity pull oras://ghcr.io/MRChemSoft/mrchem/mrchem_nompi:latest

Tagged version (here MRChem-v1.0.2 using OpenMPI-v4.0):

    $ singularity pull oras://ghcr.io/MRChemSoft/mrchem/mrchem_openmpi4.0:v1.0.2

Note that the MPI image requires that a compatible MPI library is installed and
available on the host. For information on how to launch the container:

    $ singularity run-help mrchem-mpi.sif

