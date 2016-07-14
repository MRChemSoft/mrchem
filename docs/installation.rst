============
Installation
============


-------------------
Build prerequisites
-------------------

On Stallo the supplied setup script should be able to configure things
correctly, provided all the necessary modules have been loaded::

    $ module load eigen cmake mkl boost python


-------------------------------
Obtaining and building the code
-------------------------------

An experimental version of MRChem is available on GitHub. There is no official
MRChem program yet, but you can write your own main program in a file called
``mrchem.cpp`` in the ``/pilot`` directory. You will find a sample code called
``mrchem.cpp.sample`` in this directory where some of the functionality is
demonstrated. To activate it, rename it ``mrchem.cpp`` *before* you run the
setup script::

    $ git clone git@github.com:MRChemSoft/mrchem.git
    $ cd mrchem/pilot
    $ cp mrchem.cpp.sample mrchem.cpp
    $ cd ..
    $ ./setup --release --enable-mpi=off

The setup script will now create a ``/build/pilot`` directory where you can
build the code::

    $ cd build/pilot
    $ make install -j

Feel free to do whatever you like with your pilot code, but it is your own
personal playground so don't add this file to git.


-------------------
Running the program
-------------------

By following the instructions above the code will be compiled in both MPI and
OpenMP parallel. At the moment it is only recommended to run in OpenMP parallel,
and you should use as many threads as you can spare (``unset OMP_NUM_THREADS``
will give you all threads available, otherwise use
``export OMP_NUM_THREADS N``)::

    $ unset OMP_NUM_THREADS
    $ ./mrchem-pilot.x

At the moment, no input parsing is available.
