"""
HPCCM recipe for MRChem Singularity image (MPI+OpenMP)

Contents:
  Ubuntu 20.04
  GNU compilers (upstream)
  OFED/MOFED
  PMI2 (SLURM)
  UCX
  OpenMPI 4.0.5
  CMake 3.20.6
  Python3.9
  MRChem (current source version)

Generating recipe (stdout):
  $ hpccm --recipe recipe_openmpi4.0.py --format singularity --singularity-version=3.2
"""

os_version="20.04"
cmake_version="3.20.6"
openmpi_version="4.0.5"

# CentOS base image
Stage0 += baseimage(image=f"ubuntu:{os_version}", _as='build')

# GNU compilers
compiler = gnu()
Stage0 += compiler
Stage0 += packages(
    apt=[
        "git",
        "apt-transport-https",
        "ca-certificates",
    ]
)

# (M)OFED
Stage0 += mlnx_ofed()

# UCX
Stage0 += ucx(cuda=False, ofed=True)

# PMI2
Stage0 += slurm_pmi2()

# OpenMPI (use UCX instead of IB directly)
Stage0 += openmpi(cuda=False,
                  infiniband=False,
                  pmi="/usr/local/slurm-pmi2",
                  ucx="/usr/local/ucx",
                  toolchain=compiler.toolchain,
                  version=openmpi_version)

# CMake
Stage0 += cmake(eula=True, version=cmake_version)

# Python 3
Stage0 += python(python2=False, python3=True)

# MRChem
Stage0 += raw(singularity=
"""
# copy the source tree into the container
%files
    . /mrchem

%post
    update-ca-certificates
    cd /mrchem && mkdir build
    cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/usr/local/mrchem -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpicxx -DENABLE_MPI=ON -DENABLE_OPENMP=ON -DENABLE_ARCH_FLAGS=OFF
    cmake --build build --target all -- -j$(nproc)
    cmake --build build --target install -- -j$(nproc)
"""
)


# Runtime distributable stage
Stage1 += baseimage(image=f"ubuntu:{os_version}")
Stage1 += Stage0.runtime()
Stage1 += raw(singularity=
"""
%files from build
    /usr/local/mrchem /usr/local/mrchem
"""
)

Stage1 += environment(variables={"PATH": "$PATH:/usr/local/mrchem/bin"})
Stage1 += runscript(commands=["mrchem"])

Stage1 += label(metadata={
    "Author": "Stig Rune Jensen <stig.r.jensen@uit.no>",
    "Description": "MRChem program (MPI+OpenMP version)",
    "Dependency": f"OpenMPI v{openmpi_version}"
})

help_str="""
%help
    Hybrid parallel (MPI + OpenMP) build of MRChem using OpenMPI-{openmpi_version} on a
    Ubuntu-{os_version} base image. Requires compatible OpenMPI version on the host.
    The image includes Mellanox OFED, UCX and PMI2 for compatibility with
    common HPC environments with InfiniBand and SLURM.

    For a pure OpenMP run (n threads on one process) you can run the container
    just as the regular mrchem executable, here with input file molecule.inp:

        $ export OMP_NUM_THREADS=n
        $ ./<image-name>.sif molecule

    In order to run with more that one MPI process you must first manually run
    the input parser to obtain the JSON input file. This is done by dry-running
    (--dryrun) the container on the main input file, here called molecule.inp:

        $ ./<image-name>.sif --dryrun molecule

    This will produce a new file molecule.json in the current directory which can
    be passed to the mrchem.x program inside the container using the singularity
    exec command:

        $ singularity exec <image-name>.sif mrchem.x mrchem.json

    To run in hybrid parallel (n threads on N processes) you should launch the
    singularity execution with mpirun/srun:

        $ export OMP_NUM_THREADS=n
        $ mpirun -np N singularity exec <image-name>.sif mrchem.x mrchem.json
"""
Stage1 += raw(singularity=help_str)
