rm -rf ~/Work/mrchem_forked/build/*
#./setup --prefix=~/Work/mrchem_install --omp --cxx=mpicxx --chemtensor --cmake-options="-DCMAKE_PREFIX_PATH=/home/ge84jaj/Work/mrchem_install -DCMAKE_CXX_FLAGS=-I/home/ge84jaj/Work/mrchem_install/include/MRCPP" build # --mpi
./setup --prefix=~/Work/mrchem_install --omp --mpi --cxx=mpicxx --chemtensor --cmake-options="-DCMAKE_DISABLE_FIND_PACKAGE_MRCPP=ON" build 
#./setup --prefix=~/Work/mrchem_install --omp --mpi --cxx=mpicxx --chemtensor build 

cd build
cmake ../
make -j
make install
