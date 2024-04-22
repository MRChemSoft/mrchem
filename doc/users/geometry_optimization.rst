-------------------------------
Running a geometry optimization
-------------------------------
In the following we will assume to have a valid user input file for the water
molecule called ``h2o.inp``, e.g. like this

.. include:: h2o_geopt.inp
    :literal:

A geometry optimization can be run by adding ``GeometryOptimizer`` section to any normal ``.inp`` file and setting the run keyword to ``true``:

.. include:: geopt_section.inp
    :literal:

This will start a geometry optimization with the default settings.

Obtaining accurate forces
-------------------------

In the above H_2O input example the ``world_prec`` parameter is chosen really small. This is necessary to get accurate forces. 
If a looser precision it is chosen, the geometry optimization may not converge. Pay attention to the warning::

    WARNING: Noise in force is larger than 0.2 times the larges force component!!!

    Geometry optimization onvergence cannot be guaranteed!!!

This is printed when the noise level is too high. Usually, geometry optimizations will not converge when this warning is printed.
In that case, either tighten the ``world_prec``, ``orb_thrs`` (or both) _or_ loosen the convergence criterion of the geometry optimization.


Pre-relax input geometries
--------------------------

Running high precision multi resolution wavelet calculations is computationally expensive.
It is therefore not advisable to use an input geometry with high forces, a small ``world_prec`` and start the simulation.
An optimized workflow would look something like this:

1. Optimize the geometry with a gaussian basis set. This can be done with a number of gaussian basis set codes
2. Use inaccurate forces (``world_prec`` ~ 1e-4) and a rather loose convergence criterion (``max_force_component`` ~ 1e-2) for the
   geometry optimization for a pre-relaxation with MRChem.
3. Do a tight geometry optimization (``max_force_component`` ~ 5e-4) and with an accurate MRChem calculation (``world_prec`` ~ 1e-6)


Reuse orbitals
--------------

For tight geometry optimizations where the input structure is already close to the local minimum (using cheaper pre-relaxations), it makes sense
to use the orbitals from the geometry optimization iteration *i* for the start of iteration *i+1*. This feature can be enabled by setting::

    use_previous_guess = true


Choosing an initial step size
-----------------------------

If there are some problems in the first couple of geometry optimization iterations (energy and force norm increasing) the initial step size should be chosen
manually. If a conservative choice (``init_step_size`` ~ 0.8 ) does not solve the problem, the problem is usually in the input 
geometry (wrong units, unphysical, ...) or in the potential energy surface (too much noise, error in the DFT input section, ...).

Convergence problem can be analyzed by visualizing the optimization trajectory and plots of the energy and force norm versus the geometry optimization iterations.
