Testing
-------

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
