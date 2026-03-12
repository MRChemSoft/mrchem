---------------------------------
Pseudopotential calculations
---------------------------------

MRChem supports Goedecker-Teter-Hutter (GTH) pseudopotentials, which replace
core electrons with an effective potential. This reduces the number of electrons
that need to be treated explicitly, which can significantly lower the
computational cost..

Setting up a pseudopotential calculation
----------------------------------------

A pseudopotential calculation is set up by adding a ``Pseudopotential`` section
to the input file. The ``pp_files`` keyword is a JSON string that maps elements
(or atom indices) to pseudopotential parameter files. Here is an example for a
sodium atom:

.. include:: na_pp.inp
    :literal:

The pseudopotential parameter file (here ``psppar.Na``) must be present in
the working directory when the calculation is launched. These files follow the
PSPPAR format used by BigDFT and contain the local potential
parameters, non-local projectors and optionally non-linear core correction
(NLCC) data.

Assigning pseudopotentials
--------------------------

The ``pp_files`` keyword accepts a JSON dictionary. Keys can be either element
symbols or atom indices (1-based). This allows flexible assignment:

**By element** -- all atoms of that element use the same pseudopotential::

    $pp_files
    {
    "Na": "psppar.Na",
    "Cl": "psppar.Cl"
    }
    $end

**By index** -- specific atoms can be assigned individually::

    $pp_files
    {
    "1": "psppar.Na",
    "2": "psppar.Cl"
    }
    $end

**Mixed all-electron and pseudopotential** -- setting an empty string or
``"none"`` for an element will use the all-electron description for that
element::

    $pp_files
    {
    "Na": "psppar.Na",
    "O": "none"
    }
    $end

Atoms that are not mentioned in ``pp_files`` at all will also default to
all-electron treatment.

Pseudopotential precision
-------------------------

The precision of the pseudopotential construction is controlled by
``pp_prec``, which defaults to ``world_prec``. In most cases the default
is sufficient, but it can be set explicitly if needed::

    Pseudopotential {
    $pp_files
    {
    "Na": "psppar.Na"
    }
    $end
      pp_prec = 1.0e-6
    }

PSPPAR file format
------------------

The pseudopotential parameter files use the format from the BigDFT/Goedecker
pseudopotential library. A typical file looks like this::

    spin-polarized 1       2.840       5.800  ...
     11.000  1.000  20240517                        zatom, zion, date
       12   -101130  1 0 2002 0                     pspcod, ixc, lmax, ...
        0.7499...    2 -0.7392...  -0.1720...       rloc nloc c1 .. cnloc
       2                                            nsep
       0.7125...  2  0.3650...  -0.1023...          s-projector
                         0.7522...
       0.9000...  1  0.4314...                      p-projector
                         0.0000...
        0.5941...     0.3564...                     rcore, qcore (nlcc)

The key quantities read from this file are:

- ``zion``: ionic charge (number of valence electrons)
- ``rloc``, ``nloc``, ``c1..cnloc``: local pseudopotential parameters
- ``nsep``: number of non-local (separable) projector channels
- For each channel: projector radius, dimension, and h-matrix elements
- Optional NLCC parameters: ``rcore`` and ``qcore``

GTH pseudopotential files for most elements can be obtained from the
`CP2K GTH pseudopotential repository <https://github.com/cp2k/cp2k/tree/master/data>`_.
or from here: https://github.com/OpenCPMD/GTH-pseudopotentials/tree/main
