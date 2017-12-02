============
Installation
============


-------------------
Build prerequisites
-------------------

On Stallo the supplied setup script should be able to configure things
correctly, provided all the necessary modules have been loaded. Using the Intel
tool chain (2016b does not work)::

    $ module load intel/2015b
    $ module load CMake/3.3.2-intel-2015b
    $ module load Boost/1.59.0-intel-2015b-Python-2.7.10
    $ module load Eigen/3.1.1

Using the GNU tool chain::

    $ module load foss/2016b
    $ module load CMake/3.6.2-foss-2016b
    $ module load Boost/1.61.0-foss-2016b
    $ module load Eigen/3.1.1

-------------------------------
Obtaining and building the code
-------------------------------

The public version of MRChem (currently with limited features) is available on
GitHub::

    $ git clone git@github.com:MRChemSoft/mrchem.git

To build the code using OpenMP parallelization (there is a ``--mpi`` analogue,
but this is highly experimental at the moment) with the GNU tool chain::

    $ cd mrchem
    $ ./setup --omp --prefix=.
    $ cd build
    $ make
    $ make install

With the Intel tool chain you need to specify the compilers in the setup::

    $ ./setup --cc=icc --cxx=icpc --omp --prefix=.

The official MRChem main program is located in ``src/mrchem/mrchem.cpp``, whose
executable will be built in ``build/bin/mrchem.x``. Please do not change this
file unless you know what you are doing. To try out your own ideas you can
instead write a separate main program in a file
called ``mrchem.cpp`` in the ``pilot`` directory. You will find a sample code
called ``mrchem.cpp.sample`` in this directory where some of the functionality
of MRCPP is demonstrated. To activate it, rename it ``mrchem.cpp`` *before* you
run the setup script::

    $ git clone git@github.com:MRChemSoft/mrchem.git
    $ cd mrchem/pilot
    $ cp mrchem.cpp.sample mrchem.cpp
    $ cd ..
    $ ./setup --cc=icc --cxx=icpc --omp --prefix=.
    $ cd build
    $ make
    $ make install

The pilot executable will now be built in ``build/pilot/mrchem-pilot.x``.
Feel free to do whatever you like with your pilot code, but it is your own
personal playground so don't add this file to git.

A test suite is provided (with the ``--enable-tests`` flag to setup) to make
sure that everything compiled properly::

    $ cd build/tests
    $ ./unit-tests


-------------------
Running the program
-------------------

A Python input parser will be provided along with the mrchem
executables both for the official program and the pilot code::

    $ build/bin/mrchem              // Python input parser
    $ build/bin/mrchem.x            // MRChem executable

    $ build/pilot/mrchem            // Python input parser
    $ build/pilot/mrchem-pilot.x    // Pilot executable

The input parser takes a single file argument (default ``mrchem.inp``),
processes the input and calls the main executable. Output is written to stdout
but can be redirected to an output file::

    $ ./mrchem mrchem.inp > mrchem.out &

A sample input file is provided for the pilot code. For the official program,
please refer to the MRChem manual or the examples directory.

By following the instructions above the code will be compiled in OpenMP
parallel. To run the program in parallel use the environment variable
``OMP_NUM_THREADS`` (``unset OMP_NUM_THREADS`` will give you all threads
available, otherwise use ``export OMP_NUM_THREADS N``)::

    $ export OMP_NUM_THREADS 16
    $ ./mrchem


--------
Examples
--------

There are a few examples (input ``mrchem.inp`` and expected output
``mrchem.out``) located in the ``examples`` directory. More details regarding
the input file can be found in the MRChem manual. To run e.g. the H2 molecule
(on Stallo, please **don't** do this on the login nodes)::

    $ cd examples/h2
    $ unset OMP_NUM_THREADS
    $ ../../build/bin/mrchem mrchem.inp

This will print the output from the calculation to the terminal, which you can
compare to the reference ``mrchem.out``.
