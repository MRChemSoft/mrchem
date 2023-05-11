----------------------------
Running MRChem with QCEngine
----------------------------

MRChem >=1.0 can be used as a computational engine with the `QCEngine
<http://docs.qcarchive.molssi.org/projects/qcengine/>`_ program executor.
QCEngine can be useful for running calculations on large sets of molecules and input parameters.
The results are collected in standardised `QCScheme format
<http://molssi-qc-schema.readthedocs.io/en/latest/index.html#>`_, which makes it
easy to build post-processing pipelines and store data according to Findability,
Accessibility, Interoperability, and Reuse (FAIR) of digital assets principles.
Furthermore, QCEngine provides different geometry optimization drivers that can
use the molecular gradient computed by MRChem for structural optimization.

Installation
------------

The easiest way is to install both QCEngine and MRChem in a Conda environment
using the precompiled version:

.. code-block:: bash

   conda create -n mrchem-qcng mrchem qcengine qcelemental geometric optking pip -c conda-forge
   conda activate mrchem-qcng
   python -m pip install -U pyberny

It is also possible to use your own installation of MRChem: just make sure that
the installation folder is in your ``PATH``.


.. note::

   If you want to use the precompiled, MPI-parallel version of MRChem with
   OpenMPI, install ``mrchem=*=*openmpi*`` insted of just ``mrchem``.
   A binary package compiled against MPICH is also available:
   ``mrchem=*=*mpich*``.


Single compute
~~~~~~~~~~~~~~

Calculations in QCEngine are defined in Python scripts. For example, the following runs MRChem to obtain the energy of water:

.. code-block:: python

   import qcelemental as qcel
   import qcengine as qcng


   mol = qcel.models.Molecule(geometry=[[0, 0, 0], [0, 1.5, 0], [0, 0, 1.5]],
                              symbols=["O", "H", "H"],
                              connectivity=[[0, 1, 1], [0, 2, 1]])
   print(mol)

   computation = {
       "molecule": mol,
       "driver": "energy",
       "model": {"method": "HF"},
       "keywords": {"world_prec": 1.0e-3},
   }
   ret = qcng.compute(computation, "mrchem")

   print(f"E_HF = {ret.return_result} Hartree")

You can save this sample as `mrchem-run-hf.py` and execute it with:

.. code-block:: bash

   python mrchem-run-hf.py

Which will print to screen:

.. code-block:: bash

   Molecule(name='H2O', formula='H2O', hash='b41d0c5')
   E_HF = -75.9789291596064 Hartree

Note that:

#. The molecule is specified, in Angstrom, using a `QCElemental
   <http://docs.qcarchive.molssi.org/projects/qcelemental/en/latest/>`_ object.
#. The computation is described using a Python dictionary.
#. The ``driver`` selects the kind of calculation you want to run with MRChem.
   Available drivers are: 

   - ``energy``, for single-point energy calculations.
   - ``gradient``, for evaluation of the molecular gradient at a given
     geometry.
   - ``properties``, for the calculation of molecular properties.
#. The ``model`` selects the wavefunction: ``HF`` for Hartree-Fock and any of
   the DFT functionals known to MRChem for a corresponding DFT calculation.
#. The ``keywords`` key in the dictionary accepts a dictionary of MRChem
   options. Any of the options in the usual input file are recognized.

Once you have a dictionary defining your computation, you can run it with:

.. code-block:: python

   ret = qcng.compute(computation, "mrchem")

You can reuse the same dictionary with *multiple* computational engine, *e.g.*
other quantum chemistry programs that are recognized as executors by QCEngine.
The return value from the ``compute`` function contains all data produced
during the calculation in QCSchema format including, for example, the execution
time elapsed. The full JSON output produced by MRChem is also available and can
be inspected in Python as:

.. code-block:: python

   mrchem_json_out = ret.extras["raw_output"]["output"]

The full, human-readable input is saved as the ``stdout`` property of the
object returned by ``compute``.

Parallelism
~~~~~~~~~~~

QCEngine allows you to exploit available parallel hardware.
For example, to use 20 OpenMP threads in your MRChem calculation you would
provide an additional task configuration dictionary as a ``task_config``
argument to ``compute``:

.. code-block:: python

   ret = qcng.compute(
           computation, 
           "mrchem",
           task_config={"ncores": 20}) 

You can inspect how the job was launched by printing out the ``provenance`` dictionary:

.. code-block:: python

   print(ret.extras["raw_output"]["output"]["provenance"])

.. code-block:: bash

   {
    "creator": "MRChem",
    "mpi_processes": 1,
    "routine": "/home/roberto/miniconda3/envs/mrchem-qcng/bin/mrchem.x",
    "total_cores": 1,
    "version": "1.1.0",
    "ncores": 12,
    "nnodes": 1,
    "ranks_per_node": 1,
    "cores_per_rank": 12,
    "total_ranks": 1
   }

   
It is also possible to run MPI-parallel and hybrid MPI+OpenMP jobs. Assuming
that you installed the MPICH version of the MRChem MPI-parallel Conda package,
the basic ``task_config`` argument to ``compute`` would look like:

.. code-block:: python

   task = {
     "nnodes": 1,  # number of nodes
     "ncores": 12,  # number of cores per task on each node
     "cores_per_rank": 6,  # number of cores per MPI rank
     "use_mpiexec": True,  # launch with MPI
     "mpiexec_command": "mpiexec -n {total_ranks}",  # the invocation of MPI
   }

This task configuration will launch a MPI job with 2 ranks on a single node.
Each rank has access to 6 cores for OpenMP parallelization. The ``provenance``
dictionary now shows:

.. code-block:: bash

   {
    "creator": "MRChem",
    "mpi_processes": 2,
    "routine": "mpiexec -n 2 /home/roberto/miniconda3/envs/mrchem-qcng/bin/mrchem.x",
    "total_cores": 12,
    "version": "1.1.0",
    "ncores": 12,
    "nnodes": 1,
    "ranks_per_node": 2,
    "cores_per_rank": 6,
    "total_ranks": 2
   }


The ``mpiexec_command`` is a string that will be interpolated to provide the
exact invocation. In the above example, MRChem will be run with:

.. code-block:: bash

   mpiexec -n 2 /home/roberto/miniconda3/envs/mrchem-qcng/bin/mrchem.x

The following interpolation parameters are understood by QCEngine when creating
the MPI invocation:

- ``{nnodes}``: number of nodes.
- ``{cores_per_rank}``: number of cores to use for each MPI rank.
- ``{ranks_per_node}``: number of MPI ranks per node. Computed as ``ncores // cores_per_rank``.
- ``{total_ranks}``: total number of MPI ranks. Computed as ``nnodes * ranks_per_node``.

More complex MPI invocations are possible by setting the appropriate
``mpiexec_command`` in the task configuration. For usage with a scheduler, such
as SLURM, you should refer to the documentation of your computing cluster and
the documentation of QCEngine.


Geometry optimizations
~~~~~~~~~~~~~~~~~~~~~~

Running geometry optimizations is just as easy as single compute. The following
example optimizes the structure of water using the SVWN5 functional with MW4.
The `geomeTRIC <https://geometric.readthedocs.io/en/latest/>`_ package is used
as optimization driver, but `pyberny
<https://jhrmnn.github.io/pyberny/algorithm.html>`_ or `optking
<https://optking.readthedocs.io/en/latest/?badge=latest>`_ would also work.

.. warning::

   The computation of the molecular gradient can be affected by significant
   numerical noise for MW3 and MW4, to the point that it can be impossible to
   converge a geometry optimization. Using a tighter precision might help, but
   the cost of the calculation might be prohibitively large.

.. code-block:: python

   import qcelemental as qcel
   import qcengine as qcng

   mol =  qcel.models.Molecule(
       geometry=[
           [ 0.29127930, 3.00875625, 0.20308515], 
           [-1.21253048, 1.95820900, 0.10303324],
           [ 0.10002049, 4.24958115,-1.10222079]
       ],
       symbols=["O", "H", "H"],
       fix_com=True,
       fix_orientation=True,
       fix_symmetry="c1")

   opt_input =  {
       "keywords": {
           "program": "mrchem",
           "maxiter": 70
       },
       "input_specification": {
           "driver": "gradient",
           "model": {
               "method": "SVWN5",
           },
           "keywords": {
               "world_prec": 1.0e-4,
               "SCF": {
                   "guess_type": "core_dz",
               }
           }
       },
       "initial_molecule": mol,
   }

   opt = qcng.compute_procedure(
           opt_input,
           "geometric",
           task_config={"ncores": 20})

   print(opt.stdout)

   print("==> Optimized geometry <==")
   print(opt.final_molecule.pretty_print())

   print("==> Optimized geometric parameters <==")
   for m in [[0, 1], [0, 2], [1, 0, 2]]:
       opt_val = opt.final_molecule.measure(m)
       print(f"Internal degree of freedom {m} = {opt_val:.3f}")

Running this script will print all the steps taken during the structural optimization.
The final printout contains the optimized geometry:

.. code-block:: bash

   Geometry (in Angstrom), charge = 0.0, multiplicity = 1:

      Center              X                  Y                   Z
   ------------   -----------------  -----------------  -----------------
   O                -4.146209038013     2.134923126314    -3.559202294678
   H                -4.906566693905     1.536801624016    -3.587431156799
   H                -4.270830051398     2.773072094238    -4.275607223691

and the optimized values of bond distances and bond angle:

.. code-block:: bash

   Internal degree of freedom [0, 1] = 1.829
   Internal degree of freedom [0, 2] = 1.828
   Internal degree of freedom [1, 0, 2] = 106.549
