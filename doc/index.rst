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
Sciences <https://www.mn.uio.no/hylleraas/english/>`_ at
`UiT - The Arctic University of Norway <http://en.uit.no>`_.

--------------------------------------------------------------------------------

The code is under active development, and the latest stable releases as well as
development versions can be found on `GitHub <https://github.com/MRChemSoft/mrchem>`_.

Features in MRChem-1.1:
-------------------------

* Wave functions:
    + Kohn-Sham DFT
        - Spin-polarized
        - Spin-unpolarized
        - LDA, GGA and hybrid functionals
    + Hartree-Fock
        - Restricted closed-shell
        - Unrestricted
    + Explicit external fields
        - Electric field
    + Solvent effects
        - Cavity-free PCM
* Properties:
    + Ground state energy
    + Dipole moment
    + Quadrupole moment
    + Polarizability
    + Magnetizability
    + NMR shielding constant
    + Geometric derivative
* Parallel implementation:
    + Shared memory (OpenMP): ~20 cores
    + Distributed memory (MPI): ~1000 procs
    + Hybrid scheme (MPI + OpenMP): ~10 000 cores
* Current size limitations:
    + ~2000 orbitals on ~100 high-end compute nodes (128 core/256GiB mem)
    + ~100 orbitals on a single high-memory (1TB) compute node

Upcoming features:
------------------

* Wave functions:
    + Meta-GGAs
    + ZORA Hamiltonian
    + Periodic Boundary Conditions
    + External magnetic field
* Properties:
    + Optical rotation
    + Spin-spin coupling constant
    + Hyperfine coupling constant
    + Magnetically induced currents
    + Hyperpolarizability
    + Geometry optimization
* Performance:
    + Reduced memory footprint
    + Improved DFT scaling and performance


.. toctree::
   :maxdepth: 2

   installation
   users/manual
   programmers/manual
