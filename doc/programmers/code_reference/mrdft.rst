MRDFT
=====

The classes that evaluate exchange-correlation (XC) functionals on the multiwavelet grid

Overview
--------

The XC potential is built in a few steps:

1. ``Factory`` creates the XC library backend (``XCFun`` or ``Libxc``, both
   derived from ``XCLib``) from the library name (``auto`` picks XCFun if it
   was compiled in, otherwise Libxc), and stores the functionals, the
   derivative order and the density cutoff.
2. ``Factory::build()`` asks the backend which functional families are used
   and creates the matching ``Functional`` (``LDA`` or ``GGA``, each
   spin-restricted or spin-polarized). The backend is moved into the functional, so ``build()`` can
   be called only once per ``Factory``. The functional is wrapped in an
   ``MRDFT`` object together with the ``Grid``.
3. ``MRDFT::evaluate()`` unifies the grids of the input densities and loops
   over the nodes of the grid (in parallel with OpenMP and MPI).
4. For each node ``Functional::makepot()`` transforms the density (and its
   gradient, for GGAs) to grid-point values, evaluates the XC library through
   ``XCLib::callLibEval()``, contracts the result with the higher order
   densities (``Functional::contract()``) and transforms back to the
   multiwavelet representation.

Derivative operator
-------------------

``makepot()`` differentiates the densities node by node
(``mrcpp::DerivativeCalculator::calcNode``), which only applies the
bandwidth-0 part of the derivative operator. ``ABGVOperator`` has bandwidth 0
only for a = b = 0, so ``abgv_00`` is the only valid choice. Other operators
(``bspline``, ``abgv_55``) would silently drop the terms from neighbouring
nodes. ``Factory::build()`` therefore aborts for anything but ``abgv_00``, and
the input keyword ``Derivatives.dft`` accepts only ``abgv_00``.

Adding a new XC library
-----------------------

1. Subclass ``XCLib`` (see ``XCFun`` and ``Libxc`` in
   ``src/mrdft/xc_libraries``) and implement its pure virtual functions,
   including ``name()``.
2. Add the new source file to ``src/mrdft/xc_libraries/CMakeLists.txt``,
   guarded by an ``ENABLE_<LIB>`` option like the existing libraries.
3. Add a branch for the new library name in ``Factory::Factory`` (and
   in the ``auto`` choice if appropriate), guarded by ``DISABLE_<LIB>`` as
   for the other libraries.
4. Define ``ENABLE_<LIB>`` and ``DISABLE_<LIB>`` in ``cmake/custom/main.cmake``
   (see how ``ENABLE_XCFUN`` and ``DISABLE_XCFUN`` are set), and add the
   fetch/find logic under ``external/upstream``.
5. Add the new name to the predicate of ``xc_library`` in ``python/template.yml``
   and update the generated parser and documentation.

Factory
-------
.. doxygenclass:: mrdft::Factory
   :project: MRChem
   :members:
   :protected-members:
   :private-members:

MRDFT
-----
.. doxygenclass:: mrdft::MRDFT
   :project: MRChem
   :members:
   :protected-members:
   :private-members:

Functional
----------
.. doxygenclass:: mrdft::Functional
   :project: MRChem
   :members:
   :protected-members:
   :private-members:

XCLib
-----
.. doxygenclass:: mrdft::XCLib
   :project: MRChem
   :members:
   :protected-members:
   :private-members:

XCFun
-----
.. doxygenclass:: mrdft::XCFun
   :project: MRChem
   :members:
   :protected-members:
   :private-members:

Libxc
-----
.. doxygenclass:: mrdft::Libxc
   :project: MRChem
   :members:
   :protected-members:
   :private-members:
