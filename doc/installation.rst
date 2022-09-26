============
Installation
============

-------------------
Build prerequisites
-------------------

- Python-3.7 (or later)
- CMake-3.14 (or later)
- GNU-5.4 or Intel-17 (or later) compilers (C++14 standard)

.. hint::
    We have collected the recommended modules for the different Norwegian HPC
    systems under ``tools/<machine>.env``. These files can be sourced in order
    to get a working environment on the respective machines, and may also serve
    as a guide for other HPC systems.


C++ dependencies
----------------

The MRChem program depends on the following C++ libraries:

- Input handling: `nlohmann/json-3.6  <https://github.com/nlohmann/json>`_
- Multiwavelets: `MRCPP-1.4  <https://github.com/MRChemSoft/mrcpp>`_
- Linear algebra: `Eigen-3.4  <https://gitlab.com/libeigen/eigen>`_
- DFT functionals: `XCFun-2.0  <https://github.com/dftlibs/xcfun>`_

All these dependencies will be downloaded automatically at configure time by
CMake, but can also be linked manually by setting the variables::

    MRCPP_DIR=<path_to_mrcpp>/share/cmake/MRCPP
    XCFun_DIR=<path_to_xcfun>/share/cmake/XCFun
    Eigen3_DIR=<path_to_eigen3>/share/eigen3/cmake
    nlohmann_json_DIR=<path_to_nlohmann_json>


Python dependencies
-------------------

**Users** only need a Python3 interpreter, which is used for configuration
(``setup`` script) as well as launching the program (``mrchem`` script).

**Developers** will need some extra Python packages to update the input
parser and build the documentation locally with Sphinx.

We **strongly** suggest not to install these Python dependencies globally, but
rather to use a local virtual environment. We provide a ``Pipfile`` for
specifying the Python dependencies.
We recommend using `Pipenv <https://pipenv.readthedocs.io/en/latest/>`_, since
it manages virtual environment and package installation seamlessly.
After installing it with your package manager, run::

    $ pipenv install --dev

to create a virtual environment with all developer packages installed.

The environment can be activated with::

    $ pipenv shell

Alternatively, any Python command can be run within the virtual environment by
doing::

    $ pipenv run python -c "print('Hello, world')"


-------------------------------
Obtaining and building the code
-------------------------------

The latest development version of MRChem can be found on the ``master``
branch on GitHub::

    $ git clone https://github.com/MRChemSoft/mrchem.git

The released versions can be found from Git tags ``vX.Y.Z`` under the
``release/X.Y`` branches in the same repository, or a zip file can be
downloaded from `Zenodo <https://doi.org/10.5281/zenodo.3606658>`_.

By default, all dependencies will be **fetched** at configure time if they are
not already available.


Configure
---------

The ``setup`` script will create a directory called ``<build-dir>`` and run
CMake. There are several options available for the setup, the most
important being:

``--cxx=<CXX>``
  C++ compiler [default: g++]
``--omp``
  Enable OpenMP parallelization [default: False]
``--mpi``
  Enable MPI parallelization [default: False]
``--type=<TYPE>``
  Set the CMake build type (debug, release, relwithdebinfo, minsizerel) [default: release]
``--prefix=<PATH>``
  Set the install path for make install [default: '/usr/local']
``--cmake-options=<STRING>``
  Define options to CMake [default: '']
``-h --help``
  List all options

The code can be built with four levels of parallelization:

 - no parallelization
 - only shared memory (OpenMP)
 - only distributed memory (MPI)
 - hybrid OpenMP + MPI

.. note::
    In practice we recommend the **shared memory version** for running on your
    personal laptop/workstation, and the **hybrid version** for running on a
    HPC cluster. The serial and pure MPI versions are only useful for debugging.

The default build is *without* parallelization and using GNU compilers::

    $ ./setup --prefix=<install-dir> <build-dir>

To use Intel compilers you need to specify the ``--cxx`` option::

    $ ./setup --prefix=<install-dir> --cxx=icpc <build-dir>

To build the code with shared memory (OpenMP) parallelization,
add the ``--omp`` option::

    $ ./setup --prefix=<install-dir> --omp <build-dir>

To build the code with distributed memory (MPI) parallelization, add the
``--mpi`` option *and* change to the respective MPI compilers (``--cxx=mpicxx``
for GNU and ``--cxx=mpiicpc`` for Intel)::

    $ ./setup --prefix=<install-dir> --omp --mpi --cxx=mpicxx <build-dir>

When dependencies are fetched at configuration time, they will be downloaded
into ``<build-dir>/_deps``. For the example of MRCPP, sources are saved into
the folders ``<build-dir>/_deps/mrcpp_sources-src`` and built into
``<build-dir>/_deps/mrcpp_sources-build``.

.. note::
    If you compile the MRCPP library manually as a separate project, the level
    of parallelization **must be the same** for MRCPP and MRChem. Similar
    options apply for the MRCPP setup, see
    `mrcpp.readthedocs.io <https://mrcpp.readthedocs.io/en/latest/>`_.


Build
-----

If the CMake configuration is successful, the code is compiled with::

    $ cd <build-dir>
    $ make


Test
----

A test suite is provided to make sure that everything compiled properly.
To run a collection of small unit tests::

    $ cd <build-dir>
    $ ctest -L unit

To run a couple of more involved integration tests::

    $ cd <build-dir>
    $ ctest -L integration


Install
-------

After the build has been verified with the test suite, it can be installed with
the following command::

    $ cd <build-dir>
    $ make install

This will install *two* executables under the ``<install-path>``::

    <install-path>/bin/mrchem       # Python input parser and launcher
    <install-path>/bin/mrchem.x     # MRChem executable

Please refer to the :ref:`User's Manual` for instructions for how to run the program.

.. hint::
    We have collected scripts for configure and build of the hybrid OpenMP + MPI
    version on the different Norwegian HPC systems under ``tools/<machine>.sh``.
    These scripts will build the current version under ``build-${version}``,
    run the unit tests and install under ``install-${version}``, e.g. to build
    version v1.0.0 on Fram::

        $ cd mrchem
        $ git checkout v1.0.0
        $ tools/fram.sh

    The configure step requires internet access, so the scripts must be run on
    the login nodes, and it will run on a single core, so it might take some
    minutes to complete. The scripts will *not* install the :ref:`Python
    dependencies`, so this must be done manually in order to run the code.

