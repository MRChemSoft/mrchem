-------------------
Running the program
-------------------

In the following we will assume to have a valid user input file for the water
molecule called ``h2o.inp``, e.g. like this

.. literalinclude:: h2o_getkw.inp

To run the calculation, pass the file name (without extension) as argument
to the ``mrchem`` script (make sure you understand the difference between the
``.inp``, ``.json`` and ``.out`` file, as described in the previous section)::

    $ mrchem h2o

This will under the hood actually do the following two steps::

    $ mrchem h2o.inp > h2o.json
    $ mrchem.x h2o.json > h2o.out

The first step includes input validation, which means that everything
that passes this step is a well-formed computation.


Dry-running the input parser
----------------------------

The execution of the two steps above can be done separately by dry-running the
parser script::

    $ mrchem --dryrun h2o

This will run only the input validation part and generate the ``h2o.json``
program input, but it will *not* launch the main executable ``mrchem.x``.
This can then be done manually in a subsequent step by calling::

    $ mrchem.x h2o.json

This separation can be useful for instance for developers or advanced users
who want to change some automatically generated input values before launching
the actual program, see :ref:`Input schema`.

Printing to standard output
---------------------------

By default the program will write to the text output file (``.out`` extension),
but if you rather would like it printed in the terminal you can add the
``--stdout`` option (then no text output file is created)::

    $ mrchem --stdout h2o

Reproducing old calculations
----------------------------

The JSON in/out file acts as a full record of the calculation, and can be
used to reproduce old results. Simply pass the JSON file once more
to ``mrchem.x``, and the ``"output"`` section will be overwritten::

    $ mrchem.x h2o.json

User input in JSON format
-------------------------

The user input file can be written in JSON format instead of the standard
syntax which is described in detail below. This is very convenient if you
have for instance a Python script to generate input files. The water
example above in JSON format reads (the ``coords`` string is not very elegant,
but unfortunately that's just how JSON works...):

.. literalinclude:: h2o_json.inp

which can be passed to the input parser with the ``--json`` option::

    $ mrchem --json h2o

.. note::

    A *user input file* in JSON format must **NOT** be confused with the JSON
    in/out file for the ``mrchem.x`` program. The file should still have a
    ``.inp`` extension, and contain all the same keywords which have to be
    validated and translated by the ``mrchem`` script into the ``.json``
    *program input file*.


Parallel execution
------------------

The MRChem program comes with support for both shared memory and distributed
memory parallelization, as well as a hybrid combination of the two. In order
to activate these capabilities, the code needs to be compiled with OpenMP
and/or MPI support (``--omp`` and/or ``--mpi`` options to the CMake ``setup``
script, see :ref:`Installation` instructions).


Shared memory OpenMP
++++++++++++++++++++

For the shared memory part, the program will automatically pick up the number
of threads from the environment variable ``OMP_NUM_THREADS``. If this variable
is *not* set it will usually default to the maximum available. So, to run
the code on 16 threads (all sharing the same physical memory space)::

    $ OMP_NUM_THREADS=16 mrchem h2o


Distributed memory MPI
++++++++++++++++++++++

In order to run a program in an MPI parallel fashion, it must be executed with an MPI
launcher like ``mpirun``, ``mpiexec``, ``srun``, etc. Note that it is only
the main executable ``mrchem.x`` that should be launched in parallel, **not**
the ``mrchem`` input parser script. This can be achieved *either* by running
these separately in a dry-run (here two MPI processes)::

    $ mrchem --dryrun h2o
    $ mpirun -np 2 mrchem.x h2o.json

*or* in a single command by passing the launcher string as argument to the
parser::

    $ mrchem --launcher="mpirun -np 2" h2o

This string can contain any argument you would normally pass to ``mpirun``
as it will be literally prepended to the ``mrchem.x`` command when the
``mrchem`` script executes the main program.


.. hint::

    For best performance, it is recommended to use shared memory *within*
    each `NUMA <https://en.wikipedia.org/wiki/Non-uniform_memory_access>`_
    domain (usually one per socket) of your CPU, and MPI across NUMA domains and
    ultimately machines. Ideally, the number of OpenMP threads should be
    between 8-20. E.g. on hardware with two sockets of 16 cores each, use
    OMP_NUM_THREADS=16 and scale the number of MPI processes by the size
    of the molecule, typically one process per ~5 orbitals or so (and
    definitely not *more* than one process per orbital).


Parallel pitfalls
-----------------

.. warning::

    Parallel program execution is not a black box procedure, and the behavior and
    efficiency of the run depends on several factors, like hardware configuration,
    operating system, compiler type and flags, libraries for OpenMP and MPI, type
    of queing system on a shared cluster, etc. Please make sure that the program
    runs correctly on *your* system and is able to utilize the computational
    resources before commencing production calculations.

Typical pitfalls for OpenMP
+++++++++++++++++++++++++++

- Not compiling with correct OpenMP support.
- Not setting number of threads correctly.
- **Hyper-threads:** the round-robin thread distribution might fill all
  hyper-threads on each core before moving on to the next physical core.
  In general we discourage the use of hyper-threads, and recommend a single
  thread per physical core.
- **Thread binding:** all threads may be bound to the same core, which means you
  can have e.g. 16 threads competing for the limited resources available on
  this single core (typically two hyper-threads) while all other cores are
  left idle.

Typical pitfalls for MPI
++++++++++++++++++++++++

- Not compiling with the correct MPI support.
- Default launcher options might not give correct behavior.
- **Process binding:** if a process is bound to a core, then all its spawned
  threads will also be bound to the same core. In general we recommend binding
  to socket/NUMA.
- **Process distribution:** in a multinode setup, all MPI processes might land
  on the same machine, or the round-robin procedure might count each core as
  a separate machine.

How to verify a parallel MRChem run
+++++++++++++++++++++++++++++++++++

- In the printed output, verify that MRCPP has actually been compiled with
  correct support for MPI and/or OpenMP::

    ----------------------------------------------------------------------

    MRCPP version         : 1.2.0
    Git branch            : master
    Git commit hash       : 686037cb78be601ac58b
    Git commit author     : Stig Rune Jensen
    Git commit date       : Wed Apr 8 11:35:00 2020 +0200

    Linear algebra        : EIGEN v3.3.7
    Parallelization       : MPI/OpenMP

    ----------------------------------------------------------------------

- In the printed output, verify that the correct number of processes and
  threads has been detected::

    ----------------------------------------------------------------------

     MPI processes         :      (no bank)                             2
     OpenMP threads        :                                           16
     Total cores           :                                           32

    ----------------------------------------------------------------------

- Monitor your run with ``top`` to see that you got the expected number of
  ``mrchem.x`` processes (MPI), and that they actually run at the expected
  CPU percentage (OpenMP)::

    PID   USER      PR  NI    VIRT    RES    SHR S   %CPU  %MEM     TIME+ COMMAND
    9502  stig      25   5  489456 162064   6628 R 1595,3   2,0   0:14.50 mrchem.x
    9503  stig      25   5  489596 162456   6796 R 1591,7   2,0   0:14.33 mrchem.x

- Monitor your run with ``htop`` to see which core/hyper-thread is being used
  by each process. This is very useful to get the correct binding/pinning of
  processes and threads. In general you want one threads per core, which means
  that every other hyper-thread should remain idle. In a hybrid MPI/OpenMP
  setup it is rather common that each MPI process becomes bound to a single
  core, which means that all threads spawned by this process will occupy the
  same core (possibly two hyper-threads). This is then easily detected with
  ``htop``.

- Perform dummy executions of your parallel launcher (``mpirun``, ``srun``, etc)
  to check whether it picks up the correct parameters from the resource manager
  on your cluster (SLURM, Torque, etc). You can then for instance report
  bindings and host name for each process::

    $ mpirun --print-rank-map hostname

  Play with the launcher options until you get it right. Note that Intel and
  OpenMPI have slightly different options for their ``mpirun`` and usually
  different behavior. Beware that the behavior can also change when you move
  from single- to multinode execution, so it is in general not sufficient to
  verify you runs on a single machine.

- Perform a small scaling test on e.g. 1, 2, 4 processes and/or 1, 2, 4 threads
  and verify that the total computation time is reduced as expected (don't
  expect 100% efficiency at any step).
