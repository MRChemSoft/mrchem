#!/bin/bash -l

if [ `hostname | grep -i login | wc -l` == 0 ]; then
   echo "This script MUST be run on the LOGIN node as it requires internet access"
   exit 1
fi

mrchem_dir="$(pwd)"
source ${mrchem_dir}/tools/betzy.env

cd ${mrchem_dir}
version=`cat VERSION`
build_dir=${mrchem_dir}/build-${version}
install_dir=${mrchem_dir}/install-${version}

if [ -d "${build_dir}" ]; then
    echo "Build directory already exsits, please remove"
    exit 1
else
    ./setup --prefix=${install_dir} --omp --mpi --cxx=mpicxx ${build_dir} && \
    cd ${build_dir} && \
    make && \
    OMP_NUM_THREADS=1 ctest -L unit --output-on-failure && \
    make install
fi

exit 0
