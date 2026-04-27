Testing
-------

Unit tests
++++++++++

We perform unit testing of our code. The unit testing framework used is
`Catch <https://github.com/philsquared/Catch>`_ The framework provides quite an
extensive set of macros to test various data types, it also provides facilities
for easily setting up test fixtures.  Usage is extremely simple and the
`documentation <https://github.com/philsquared/Catch/blob/master/docs/Readme.md>`_
is very well written.  For a quick primer on how to use Catch refer to:
https://github.com/philsquared/Catch/blob/master/docs/tutorial.md
The basic idea of unit testing is to test each building block of the code
separataly. In our case, the term "building block" is used to mean a class.

To add new tests for your class you have to:

#. create a new subdirectory inside tests/ and add a line like the following
   to the CMakeLists.txt

   .. code-block:: cmake

      add_subdirectory(new_subdir)

#. create a CMakeLists.txt inside your new subdirectory.
   This CMakeLists.txt adds the source for a given unit test to the global ``UnitTestsSources``
   property and notifies CTest that a test with given name is part of the test suite.
   The generation of the CMakeLists.txt can be managed by ``make_cmake_files.py`` Python script.
   This will take care of also setting up CTest labels. This helps in further grouping
   the tests for our convenience.
   Catch uses tags to index tests and tags are surrounded by square brackets. The Python script
   inspects the sources and extracts labels from Catch tags.
   The ``add_Catch_test`` CMake macro takes care of the rest.

   We require that each source file containing tests follows the naming convention
   new_subdir_testname and that testname gives some clue to what is being tested.
   Depending on the execution of tests in a different subdirectory is bad practice.
   A possible workaround is to add some kind of input file and create a text fixture
   that sets up the test environment. Have a look in the ``tests/input`` directory
   for an example

Running Different Types of Tests
++++++++++++++++++++++++++++++++

The tests in MRChem should all have a list of labels in the CMakeLists.txt inside their respective subdirectories, all in lower case letters. 
CTest allows you to call all tests containing a specific label by running the following command in you build directory:

   .. code-block:: bash

      ctest -L "label"    # label = scf, dft, H2, polarization, ...

You can also use the -R command to run specific tests by name or parts of a name. 
Note that command is utilizing regex, so if multiple test names contain the entered keyword, they will all run.
Also note that this is command is case sensitive, and the casing of the test name may not correspond with that of the test subdirectory.

   .. code-block:: bash

      ctest -R "name"    # name = H2O_energy_BLYP, ...


Division of Tests
+++++++++++++++++

The tests are divided into three sets, depending on the runtime cost: short, medium and long.
Which set a test belongs to is included in the labels of the test as one of the three mentioned keywords.

- Short tests are smaller critical tests (up to ~1 min.) that are run every time a new commit is pushed to github
- Medium tests are longer, and it not part of the default test pipeline in github. However, they should always be run before a pull request is merged
- Long tests are more extensive tests that should only be run sometimes

