#!/bin/bash -l

if [ -z "$1" ]; then
    # If not, prompt the user for input
    echo "This script requires an NRIS account for compilation on a compute node"
    read -p "Please state your account (nnxxxxk): " account
else
    # Otherwise, use the argument
    account="$1"
fi

mrchem_dir="$(pwd)"
source ${mrchem_dir}/tools/olivia.env

cd ${mrchem_dir}
version=`cat VERSION`
build_dir=${mrchem_dir}/build-${version}
install_dir=${mrchem_dir}/install-${version}

if [ -d "${build_dir}" ]; then
    echo "Build directory already exists, please remove"
    exit 1
else
    ./setup --prefix=${install_dir} --omp --mpi ${build_dir} && \
    cd ${build_dir} && \
    srun --cpus-per-task=4 --mem-per-cpu=2G --time=01:00:00 --account=$account bash -c \
    "source ${mrchem_dir}/tools/olivia.env ; \
    make -j4 ; \
    OMP_NUM_THREADS=4 ctest -L unit --output-on-failure ; \
    make install -j4"
fi

exit 0
