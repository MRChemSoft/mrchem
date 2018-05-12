===================
MRChem: User manual
===================

---------------------
The mrchem input file
---------------------

The input file is organized in sections and keywords that can be of different
type

.. code-block:: bash

     Section {
        keyword_1 = 1                       # int
        keyword_2 = 3.14                    # double
        keyword_3 = [1, 2, 3]               # int array
        keyword_4 = foo                     # string
        keyword_5 = true                    # bool
     }

Top section
-----------

The main input section contain the single most important keyword:
the relative precision :math:`\epsilon_{rel}` that will be guaranteed in the
calculation. The top section is not specified by name, just write the keywords
directly, e.g

.. code-block:: bash

    rel_prec = 1.0e-5                       # Overall relative precision

The relative precision sets an upper limit for the number of correct digits
you are expected to get out of the computation (note that
:math:`\epsilon_{rel}=10^{-6}` yields :math:`\mu` Ha accuracy for the hydrogen
molecule, but only mHa accuracy for benzene). It is also possible to specify
an `absolute` precision for the molecular energy by replacing ``rel_prec``
with ``abs_prec``. This will provide e.g. mHa precision in the total energy
regardless of the molecular size (this might get `very` expensive for large
systems). In this case the magnitude of the energy is estimated as

.. math:: \tilde{E} = \sum_i^{nuc} Z_i^{5/2} 

and the relative precision is set as

.. math:: \epsilon_{rel} = \frac{\epsilon_{abs}}{\tilde{E}}

With the ``abs_prec`` keyword one can also use kcal/mol or kJ/mol as energy
unit instead of Hartree by setting the ``energy_unit`` keyword. The following
will set ``rel_prec`` sufficiently high so that the energy can be computed
within a kcal/mol

.. code-block:: bash

    abs_prec    = 1.0
    energy_unit = kcal

Note that ``rel_prec`` is the fundamental input parameter that is finally
passed to the program, so if this is set explicitly in the input file, it will
always have precedence.

MRChem uses a smoothed nuclear potential to avoid numerical problems in
connection with the :math:`Z/|r-R|` singularity. The smoothing is controlled by
a single parameter ``nuc_prec`` that is related to the expected error in the
energy due to the smoothing.

.. code-block:: bash

    nuc_prec = 1.0e-6

If the ``nuc_prec`` keyword is left out it is set equal to ``rel_prec``.

A final keyword in the main input section that is usually used for debugging is
the ``printlevel`` (should be zero for production calculations).

MRA
-----

The MultiResolution Analysis (MRA) input section defines the polynomial basis
as well as the computational domain of the calculation (defaults shown except
for ``order``, see below)

.. code-block:: bash

    MRA {
        order          = 7                  # Polynomial order of MW basis
        basis_type     = Interpolating      # Legendre or Interpolating
        min_scale      = 0                  # Size of each root box 2^{-n}
        max_scale      = 20                 # Maximum level of refinement 2^{-n}
        boxes          = [ 1, 1, 1 ]        # Number of root boxes
        corner         = [ 0, 0, 0]         # Translation of first root box
        gauge_origin   = [0.0, 0.0, 0.0]    # Origin used in molecular properties
        center_of_mass = false              # Use COM gauge origin
    }

The MW basis is defined by the polynomial order :math:`k`, and the type of
scaling functions (Legendre or Interpolating polynomials). Note that
increased precision requires higher polynomial order (use e.g :math:`k = 5`
for :math:`\epsilon_{rel} = 10^{-3}`, and :math:`k = 13` for
:math:`\epsilon_{rel} = 10^{-9}`, and interpolate in between). If the ``order``
keyword is left out it will be set automatically according to

.. math:: k=-1.5*log_{10}(\epsilon_{rel})

The scale and translation of the root boxes are absolute, which means that the
only way to get a symmetric world around the origin is to use two root ``boxes``
in each direction and set ``corner`` at -1 (if this does not fit well with your
molecular geometry, use a larger box or translate your molecular coordinates).
The computational world should be large enough so that the electron density
vanishes at the boundaries. The ``gauge_origin`` can also be specified (relevant
for molecular properties), or set to the molecular ``center_of_mass``. The
default computational domain displayed above corresponds to the unit cube (in
bohr). The maximum refinement level ``max_scale`` should preferably be as small
as possible for computational efficiency, but if this level is actually
encountered in the calculation, the accuracy might be affected. Note that
the total span of length scales (``max_scale`` - ``min_scale``) cannot exceed
32 (integer precision is :math:`2^{32}`).

Molecule
--------

This input section specifies the geometry, charge and spin multiplicity of the
molecule, e.g. for water (coords must be specified, otherwise
defaults are shown)

.. code-block:: bash

    Molecule {
        charge       = 0                    # total charge of molecule
        multiplicity = 1                    # spin multiplicity
        angstrom     = false                # geometry given in angstrom
        $coords
        O   0.0000     0.0000     0.0000
        H   0.0000     1.4375     1.1500
        H   0.0000    -1.4375     1.1500
        $end
    }

WaveFunction
------------

Here we give the wavefunction method and whether we run spin restricted (alpha
and beta spins are forced to occupy the same spatial orbitals) or not (method
must be specified, otherwise defaults are shown) 

.. code-block:: bash

    WaveFunction {
        method     = <wavefunction_method>  # Core, Hartree, HF or DFT
        restricted = true                   # Spin restricted/unrestricted
    }

There are currently four methods available: Core Hamiltonian, Hartree,
Hartree-Fock (HF) and Density Functional Theory (DFT). When running DFT the
functional(s) must be specified in a separate DFT section (see below).

DFT
---
 
This section specifies the exchange-correlation functional used in DFT
(functional names must be specified, otherwise defaults are shown)

.. code-block:: bash

    DFT {
        spin = false                        # Use spin-polarized functionals
        exact_exchange = 0.0                # Amount of exact HF exchange
        density_cutoff = 0.0                # Cutoff to set XC potential to zero
        explicit_der = true                 # Use explicit derivatives or gamma
        $functionals
        <func1>     1.0                     # Functional name and coefficient
        <func2>     1.0
        $end
    }

You can specify as many functionals as you want, and they will be added on top
of each other with the given coefficient. Both exchange and correlation
functionals must be set explicitly, e.g. ``SLATERX`` and ``VWN5C`` for the
standard LDA functional. If the ``spin`` parameter is not explicitly set it will
follow the ``restricted`` parameter of the ``WaveFunction`` section.
Option to use explicit partial derivatives for the density gradients
(:math:`\delta f_{xc}/\delta\nabla\rho`) or the invariants
(:math:`\gamma=\nabla\rho\cdot\nabla\rho`). For hybrid functionals you must
specify the amount of exact Hartree-Fock exchange that should be used (0.2 for
B3LYP and 0.25 for PBE0 etc.). Option to use spin-polarized functionals.
XC functionals are provided by the `XCFun <https://github.com/dftlibs/xcfun>`_
library.

Properties
----------

Specify which properties to compute. Currently the following are available
(defaults shown)

.. code-block:: bash

    Properties {
        scf_energy    = false               # Compute total SCF energy
        dipole_moment = false               # Compute dipole moment
    }

SCF
---

Specify the parameters for the SCF optimization of the ground state wave
function (defaults shown)

.. code-block:: bash

    SCF {
        run            = true              # Run SCF optimization
        orbital_thrs   = -1.0              # Convergence threshold orbitals
        property_thrs  = -1.0              # Convergence threshold energy
        orbital_prec   = [1.0e-4, -1.0]    # Initial and final relative precision in SCF
        kain           = 0                 # Length of KAIN iterative subspace
        rotation       = 0                 # Iterations between each localization/diagonalization
        max_iter       = -1                # Maximum number of SCF iterations
        canonical      = false             # Use canonical or localized  orbitals
        write_orbitals = false             # Write final orbitals to disk
        initial_guess  = none              # Type of inital guess (none, gto, mw)
    }

With ``run=false`` no SCF optimization is performed, and the requested molecular
properties are computed directly from the initial guess wave function.

We specify a convergence threshold both for the orbitals
(:math:`\|\Delta \phi_i \|`) and the property (:math:`\Delta E`). The default
value of -1.0 means that the threshold will not be considered in the
optimization. The property (total SCF energy) should converge quadratically in
the orbital errors, however, it will still be limited by the overall precision
``rel_prec`` in the calculation. For instance, the following will converge the
energy within nine digits, but only five of them are guaranteed to be correct

.. code-block:: bash

    rel_prec = 1.0e-5

    SCF {
        property_thrs = 1.0e-9
    }

When computing other properties than total energy, the important threshold is
that for the orbitals, which translates approximately to the relative accuracy
that you can expect for other properties. The following input should give five
digits for the dipole moment (always keep a factor of 10 between ``rel_prec``
and ``orbital_thrs`` to avoid numerical instabilities)

.. code-block:: bash

    rel_prec = 1.0e-6

    SCF {
        orbital_thrs = 1.0e-5
    }

If *both* thresholds are omitted in this section they will be
set according to the top level ``rel_prec``

.. math:: \Delta E < \frac{\epsilon_{rel}}{10}
.. math:: \|\Delta \phi_i \| < \sqrt{\frac{\epsilon_{rel}}{10}}

This should yield a final energy accurate within the chosen relative precision.
This means that in order to get for instance milli-Hartree accuracy in energy,
you need only specify the ``abs_prec`` keyword in the top level, then all
related parameters (``order``, ``rel_prec``, ``nuc_prec``, ``orbital_thrs`` and
``property_thrs``) will be adjusted so that the requested precision is reached.

The ``orbital_prec=[init,final]`` keyword controls the dynamic precision used
in the SCF iterations. To improve efficiency, the first iterations are done
with reduced precision, starting at ``init`` and gradually increased
to ``final``. The initial precision should not be set lower than
``init=1.0e-3``, and the final precision should not exceed the top level
``rel_prec``. Negative values sets them equal to ``rel_prec``. 

The ``kain`` keyword sets the size of the iterative subspace that is used
in the KAIN accelerator for the orbital optimization.

The ``rotation`` and ``canonical`` keywords says how often the Fock matrix
should be diagonalized/localized (for iterations in between, a Löwdin
orthonormalization using the overlap matrix :math:`S^{-1/2}` is used).
Option to use Foster-Boys localization or Fock matrix diagonalization in
these rotations. Note that the KAIN history is cleared every time this
rotation is employed to avoid mixing of orbitals in the history, so
``rotation=1`` effectively cancels the KAIN accelerator. The default
``rotation=0`` will localize/diagonalize the first two iterations and then
perform Löwdin orthonormalizations from that point on (this is usually the
way to go).

You also need to specify which ``initial_guess`` to use, "none" means starting
from hydrogen solutions (this requires no extra input, but is a quite poor
guess), "gto" means starting with a wave function from a converged calculation
using a small GTO basis set (basis and MO matrix input files must be provided)
and "mw" means starting from a previous MRChem calculation (compatible orbitals
must have been written to disk using the ``write_orbitals`` keyword).

Example 1
---------

The following input will compute the Hartree-Fock energy of water to
micro-Hartree precision

.. code-block:: bash

    abs_prec = 1.0e-6

    MRA {
        min_scale = -5                      # Size of each root box 2^{-n}
        boxes     = [ 2, 2, 2]              # Number of root boxes
        corner    = [-1,-1,-1]              # Translation of first root box
    }

    Molecule {
        $coords
        O   0.0000     0.0000     0.0000
        H   0.0000     1.4375     1.1500
        H   0.0000    -1.4375     1.1500
        $end
    }

    WaveFunction {
        method = HF                         # Core, Hartree, HF or DFT
    }

    Properties {
        scf_energy = true                   # Compute total energy
    }

    SCF {
        kain = 3                            # Length of KAIN iterative subspace
    }


Example 2
---------

The following input will compute the B3LYP energy (six digits) and dipole moment
(four digits) of carbon monoxide 

.. code-block:: bash

    rel_prec = 1.0e-6

    MRA {
        min_scale = -5                      # Size of each root box 2^{-n}
        boxes     = [ 2, 2, 2]              # Number of root boxes
        corner    = [-1,-1,-1]              # Translation of first root box
    }

    Molecule {
        angstrom = true
        $coords
        C   0.0000     0.0000    -0.56415
        O   0.0000     0.0000     0.56415
        $end
    }

    WaveFunction {
        method = DFT                        # Core, Hartree, HF or DFT
    }

    DFT {
        exact_exchange = 0.20               # Amount of exact HF exchange
        $functionals
        BECKEX      0.80                    # Functional name and coefficient
        LYPC        1.00
        $end
    }

    Properties {
        scf_energy = true                   # Compute total energy
        dipole_moment = true                # Compute dipole moment
    }

    SCF {
        kain          = 3                   # Length of KAIN iterative subspace
        orbital_thrs  = 1.0e-4              # Convergence threshold orbitals
        property_thrs = 1.0e-7              # Convergence threshold energy
    }
