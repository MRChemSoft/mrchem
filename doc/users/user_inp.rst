---------------
User input file
---------------

The input file is organized in sections and keywords that can be of different
type. Input keywords and sections are **case-sensitive**, while `values` are
**case-insensitive**.

.. code-block:: bash

    Section {
       keyword_1 = 1                       # int
       keyword_2 = 3.14                    # double
       keyword_3 = [1, 2, 3]               # int array
       keyword_4 = foo                     # string
       keyword_5 = true                    # bool
    }

A  complete list of available input keywords can be found in the :ref:`Input
parameters`.

Top section
-----------

The main input section contain two keywords, the relative precision
:math:`\epsilon_{rel}` that will be guaranteed in the calculation and the size
of the computational domain. The top section is not specified by name, just
write the keywords directly, e.g

.. code-block:: bash

    world_prec = 1.0e-5                       # Overall relative precision
    world_size = 5                            # World size 2^world_size

The relative precision sets an upper limit for the number of correct digits
you are expected to get out of the computation (note that
:math:`\epsilon_{rel}=10^{-6}` yields :math:`\mu` Ha accuracy for the hydrogen
molecule, but only mHa accuracy for benzene).

The computational domain is always symmetric around the origin, with `total`
size given by the ``world_size`` parameter as :math:`[2^n]^3`, e.i.
``world_size=5`` gives a domain of :math:`[-16,16]^3` atomic units (bohrs).
Make sure that the world
is large enough to allow the molecular density to reach zero on the boundary. 
The ``world_size`` parameter can be left out, in which case the size will be
estimated based on the molecular geometry.

Precisions
----------

MRChem uses a smoothed nuclear potential to avoid numerical problems in
connection with the :math:`Z/|r-R|` singularity. The smoothing is controlled by
a single parameter ``nuc_prec`` that is related to the expected error in the
energy due to the smoothing. There are also different precision parameters for
the `construction` of the Poisson and Helmholtz integral operators.

.. code-block:: bash

    Precisions {
        nuc_prec = 1.0e-6
        poisson_prec = 1.0e-6
        helmholtz_prec = 1.0e-6
    }

By default, all precision parameters follows ``world_prec``.

Printer
-------

The following are the most important Printer keywords:

.. code-block:: bash

    Printer {
        print_level = 0       # Level of detail in the written output
        print_input = true    # Print the input file at the start of the calculation
    }

Available print levels:

- ``print_level=0`` prints mainly properties
- ``print_level=1`` adds timings for individual steps
- ``print_level=2`` adds memory and timing information on ``OrbitalVector`` level
- ``print_level=3`` adds memory and timing information on ``Orbital`` level
- ``print_level>10`` adds a *lot* more output from deep within MRCPP



Basis
-----

This section defines the polynomial MultiWavelet basis

.. code-block:: bash

    Basis {
        type = Interpolating  # Legendre or Interpolating
        order = 7             # Polynomial order of MW basis
    }

The MW basis is defined by the polynomial order :math:`k`, and the type of
scaling functions (Legendre or Interpolating polynomials). Note that
increased precision requires higher polynomial order (use e.g :math:`k = 5`
for :math:`\epsilon_{rel} = 10^{-3}`, and :math:`k = 13` for
:math:`\epsilon_{rel} = 10^{-9}`, and interpolate in between). If the ``order``
keyword is left out it will be set automatically according to

.. math:: k=-1.5*log_{10}(\epsilon_{rel})

The Basis section can usually safely be omitted in the input.

Molecule
--------

This input section specifies the geometry, charge and spin multiplicity of the
molecule, e.g. for water (coords must be specified, otherwise defaults are
shown)

.. code-block:: bash

    Molecule {
        charge       = 0                    # total charge of molecule
        multiplicity = 1                    # spin multiplicity
        translate    = false                # translate center of mass to origin
        angstrom     = false                # geometry given in angstrom
    $coords
    O   0.0000     0.0000     0.0000
    H   0.0000     1.4375     1.1500
    H   0.0000    -1.4375     1.1500
    $end
    }

Since the computational domain is always cubic and symmetric around the origin
it is usually a good idea to ``translate`` the molecule to the origin.

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
Hartree-Fock (HF) and Density Functional Theory (DFT). When running DFT you can
`either` set one of the default functionals in this section (e.g. ``method =
B3LYP``), `or` you can set ``method = DFT`` and specify a "non-standard"
functional in the separate DFT section (see below). See :ref:`Input
parameters` for a list of available default functionals.

DFT
---

This section specifies the exchange-correlation functional used in DFT
(functional names must be specified, otherwise defaults are shown)
This section can be omitted if you are using a default functional, see above.

.. code-block:: bash

    DFT {
        spin = false                        # Use spin-polarized functionals
        use_gamma = true                    # Use explicit derivatives or gamma
        density_cutoff = 0.0                # Cutoff to set XC potential to zero
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
specify the amount of exact Hartree-Fock exchange as a separate functional
``EXX`` (``EXX 0.2`` for B3LYP and ``EXX 0.25`` for PBE0 etc.). Option to use
spin-polarized functionals. Unrestricted calculations will use spin-polarized
functionals by default. The XC functionals are provided by the
`XCFun <https://github.com/dftlibs/xcfun>`_ library.

Properties
----------

Specify which properties to compute. Currently the following are available
(defaults shown)

.. code-block:: bash

    Properties {
        scf_energy    = true                # Compute total SCF energy
        dipole_moment = false               # Compute dipole moment
    }

SCF
---

Specify the parameters for the SCF optimization of the ground state wave
function (defaults shown)

.. code-block:: bash

    SCF {
        kain           = 0                 # Length of KAIN iterative subspace
        max_iter       = -1                # Maximum number of SCF iterations
        rotation       = 0                 # Iterations between diagonalize/localize
        localize       = false             # Use canonical or localized  orbitals
        orbital_thrs   = -1.0              # Convergence threshold orbitals
        property_thrs  = -1.0              # Convergence threshold energy
        initial_guess  = sad_dz            # Type of inital guess (mw, gto, core, sad)
    }

We specify a convergence threshold both for the orbitals
(:math:`\|\Delta \phi_i \|`) and the property (:math:`\Delta E`). The default
value of -1.0 means that the threshold will not be considered in the
optimization. The property (total SCF energy) should converge quadratically in
the orbital errors. However, it will still be limited by the overall precision
``world_prec`` in the calculation. For instance, the following will converge the
energy within nine digits, but only five of them are guaranteed to be correct

.. code-block:: bash

    world_prec = 1.0e-5

    SCF {
        property_thrs = 1.0e-9
    }

When computing other properties than total energy, the important threshold is
that for the orbitals, which translates approximately to the relative accuracy
that you can expect for other properties. The following input should give five
digits for the dipole moment (always keep a factor of 10 between ``world_prec``
and ``orbital_thrs`` to avoid numerical instabilities)

.. code-block:: bash

    world_prec = 1.0e-6

    SCF {
        orbital_thrs = 1.0e-5
    }

If *both* thresholds are omitted in this section they will be
set according to the top level ``world_prec``

.. math:: \Delta E < \frac{\epsilon_{rel}}{10}
.. math:: \|\Delta \phi_i \| < \sqrt{\frac{\epsilon_{rel}}{10}}

This should yield a final energy accurate within the chosen relative precision.

The ``kain`` keyword sets the size of the iterative subspace that is used
in the KAIN accelerator for the orbital optimization.

The ``rotation`` and ``localize`` keywords says how often the Fock matrix
should be diagonalized/localized (for iterations in between, a Löwdin
orthonormalization using the overlap matrix :math:`S^{-1/2}` is used).
Option to use Foster-Boys localization or Fock matrix diagonalization in
these rotations. Note that the KAIN history is cleared every time this
rotation is employed to avoid mixing of orbitals in the history, so
``rotation=1`` effectively cancels the KAIN accelerator. The default
``rotation=0`` will localize/diagonalize the first two iterations and then
perform Löwdin orthonormalizations from that point on (this is usually the
way to go). See :ref:`Input parameters` for more details.

Plotter
-------

It is possible to get a 3D cube plot of the converged orbitals and density by
setting the keywords ``plot_orbital`` and ``plot_density`` in the ``SCF``
section. In addition the ``Plotter`` section must be specified in order to get
correct plotting parameters. The functions are plotted in the volume spanned by
the three vectors A, B and C, relative to the origin O. The following example
will generate a 20x20x20 plot in the volume [-4,4]^3 of the density plus
orbitals 1 and 2:

.. code-block:: bash

    Plotter {
      points = [20, 20, 20]  # number of grid points
      O = [-4.0,-4.0,-4.0]   # plot origin
      A = [8.0, 0.0, 0.0]    # boundary vector
      B = [0.0, 8.0, 0.0]    # boundary vector
      C = [0.0, 0.0, 8.0]    # boundary vector
    }

    SCF {
      plot_density = true    # plot converged density (including spin for open-shell)
      plot_orbital = [1,2]   # plot converged 1 and 2 (negative idx plots all)
    }

The generated files (e.g. ``plots/phi_1_re.cube``) can be viewed directly in a
web browser by `blob <https://github.com/densities/blob/>`_ , like this benzene
orbital:

.. image:: gfx/blob.png

Example 1
---------

The following input will compute the Hartree-Fock energy of water to
six digits precision, world size :math:`[-32,32]^3`

.. code-block:: bash

    world_prec = 1.0e-6
    world_size = 6

    Molecule {
        translate = true
    $coords
    O   0.0000     0.0000     0.0000
    H   0.0000     1.4375     1.1500
    H   0.0000    -1.4375     1.1500
    $end
    }

    WaveFunction {
        method = HF
    }

    Properties {
        scf_energy = true
    }

    SCF {
        kain = 3
    }


Example 2
---------

The following input will compute the B3LYP energy (six digits) and dipole moment
(four digits) of carbon monoxide, automatic world size

.. code-block:: bash

    world_prec = 1.0e-6

    Molecule {
        angstrom = true
    $coords
    C   0.0000     0.0000    -0.56415
    O   0.0000     0.0000     0.56415
    $end
    }

    WaveFunction {
        method = B3LYP
    }

    Properties {
        scf_energy = true
        dipole_moment = true
    }

    SCF {
        kain          = 3
        orbital_thrs  = 1.0e-4
        property_thrs = 1.0e-7
    }
