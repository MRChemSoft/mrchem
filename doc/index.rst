.. MRChem documentation master file, created by
   sphinx-quickstart on Tue Jan 26 15:03:29 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==================================
Welcome to MRChem's documentation!
==================================

MRChem is a numerical real-space code for molecular
electronic structure calculations within the self-consistent field (SCF)
approximations of quantum chemistry (Hartree-Fock and Density Functional
Theory). The code is divided in two main parts: the `MultiResolution
Computation Program Package <https://mrcpp.readthedocs.io/en/latest>`_ (MRCPP),
which is a general purpose numerical
mathematics library based on multiresolution analysis and the multiwavelet
basis which provide low-scaling algorithms as well as rigorous error control
in numerical computations, and the MultiResolution Chemistry (MRChem) program
that uses the functionalities of MRCPP for computational chemistry applications.

The code is being developed at the `Hylleraas Centre for Quantum Molecular
Sciences <http://www.ctcc.no/>`_ at
`UiT - The Arctic University of Norway <http://en.uit.no>`_.

--------------------------------------------------------------------------------

We are currently in the process of rewriting the code and making it publicly
available, and the latest version (with limited functionality) can be found on
`GitHub <https://github.com/MRChemSoft/mrchem>`_. This is **not** a stable
version, expect major changes in the future.

Features as of March 2017:
--------------------------

* Wave functions:
    + Spin-unpolarized Kohn-Sham DFT (LDA, GGA and hybrid)
    + Restricted Hartree-Fock (closed-shell only)
    + Unrestricted Hartree-Fock
* Properties:
    + Ground state energy
    + Dipole moment
* Parallel implementation:
    + Shared memory (OpenMP): ~20 cores
* Current limitations on a single high-memory compute node (768GB):
    + nano-Hartree accuracy: ~10 orbitals
    + micro-Hartree accuracy: ~50 orbitals
    + milli-Hartree accuracy: ~100 orbitals

Upcoming features:
------------------

* Wave functions:
    + Spin-polarized Kohn-Sham DFT
    + Meta-GGAs
* Properties:
    + Quadrupole moment
    + Polarizability
    + Hyperpolarizability
    + Optical rotation
    + Magnetizability
    + NMR shielding constant
    + Spin-spin coupling constant
    + Hyperfine coupling constant
    + Magnetically induced currents
* Parallel implementation:
    + Distributed memory (MPI)
    + Hybrid scheme (MPI + OpenMP)


.. toctree::
   :maxdepth: 2

   installation
   mrchem_manual
   code_reference/classes-and-functions
