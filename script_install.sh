rm -rf /work/ge84jaj/mrchem_forked/build/*
./setup --prefix=/work/ge84jaj/mrchem_install --omp --mpi --cxx=mpicxx build
cd build
cmake ../
make -j
make install
