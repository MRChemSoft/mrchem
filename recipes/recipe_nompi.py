"""
HPCCM recipe for MRChem Singularity image (OpenMP)

Contents:
  Ubuntu 20.04
  GNU compilers (upstream)
  CMake 3.20.6
  Python3.9
  MRChem (current source version)

Generating recipe (stdout):
  $ hpccm --recipe recipe_nompi.py --format singularity --singularity-version=3.2
"""

os_version="20.04"
cmake_version="3.20.6"

# Ubuntu base image
Stage0 += baseimage(image=f"ubuntu:{os_version}", _as="build")

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
    cmake -B build -S . -DCMAKE_INSTALL_PREFIX=/usr/local/mrchem -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=OFF -DENABLE_OPENMP=ON -DENABLE_ARCH_FLAGS=OFF
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
    "Description": "MRChem program (OpenMP version)"
})

help_str="""
%help
    Shared memory parallel (OpenMP) build of MRChem on a Ubuntu-{os_version} base image.

    For a pure OpenMP run (n threads on one process) you can run the container
    just as the regular mrchem executable, here with input file molecule.inp:

        $ export OMP_NUM_THREADS=n
        $ ./<image-name>.sif molecule
"""
Stage1 += raw(singularity=help_str)
