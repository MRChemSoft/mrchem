CMake usage
===========

This is a brief guide to our CMake infrastructure which is managed
*via* `Autocmake <http://autocmake.readthedocs.org/en/latest/>`_

.. warning::

   The minimum required CMake version is 2.8.8

Adding new source subdirectories and/or files
---------------------------------------------

Developers **HAVE TO** manually list the sources in a given subdirectory
of the main source directory ``src/``.

New subdirectory
................

First of all, you will have to let CMake know that a new source-containing
subdirectory has been added to the source tree. Due to the hierarchical
approach CMake is based upon, you will need to modify the ``CMakeLists.txt`` in
the ``src`` directory and create a new one in your new subdirectory.  For the
first step:

   1. if your new subdirectory contains header files, add a line like
   the following to the ``CMakeLists.txt`` file contained in the ``src`` directory:

   .. code-block:: cmake

      ${CMAKE_CURRENT_LIST_DIR}/subdir_name

   to the command setting the list of directories containing headers.  This
   sets up the list of directories where CMake will look for headers with
   definitions of classes and functions. If your directory contains Fortran
   code you can skip this step;

   2. add a line like the following to the ``CMakeLists.txt`` file contained in the
   ``src`` directory:

   .. code-block:: cmake

      add_subdirectory(subdir_name)

   This will tell CMake to go look inside ``subdir_name`` for a ``CMakeLists.txt``
   containing more sets of instructions.  It is preferable to add these new
   lines in **alphabetic order**

Inside your new subdirectory you will need to add a ``CMakeLists.txt`` file containing
the set of instructions to build your cutting edge code. This is the second step.
Run the ``make_cmake_files.py`` Python script in the ``src/`` directory:

.. code-block:: bash

   python make_cmake_files.py --libname=cavity --lang=CXX

to generate a template ``CMakeLists.txt.try`` file:

.. code-block:: cmake

   # List of headers
   list(APPEND headers_list Cavity.hpp ICavity.hpp Element.hpp GePolCavity.hpp RegisterCavityToFactory.hpp RestartCavity.hpp)

   # List of sources
   list(APPEND sources_list ICavity.cpp Element.cpp GePolCavity.cpp RestartCavity.cpp)

   add_library(cavity OBJECT ${sources_list} ${headers_list})
   set_target_properties(cavity PROPERTIES POSITION_INDEPENDENT_CODE 1 )
   set_property(GLOBAL APPEND PROPERTY PCMSolver_HEADER_DIRS ${CMAKE_CURRENT_LIST_DIR})
   # Sets install directory for all the headers in the list
   foreach(_header ${headers_list})
      install(FILES ${_header} DESTINATION include/cavity)
   endforeach()

The template might need additional editing.
Each source subdirectory is the lowest possible in the CMake
hierarchy and it contains set of instructions for:

#. exporting a list of header files (.h or .hpp) to the upper level in the
   hierarchy, possibly excluding some of them
#. define install targets for the files in this subdirectory.

All the source files are compiled into the unique static library ``libpcm.a`` and unique
dynamic library ``libpcm.so``.
This library is the one the host QM program need to link.

Searching for libraries
.......................

In general, the use of the `find_package <http://www.cmake.org/cmake/help/v3.0/command/find_package.html>`_
macro is to be preferred, as it is standardized and ensured to work on any
platform.  Use of ``find_package`` requires that the package/library you want to
use has already a module inside the CMake distribution.  If that's not the
case, you should *never* use the following construct for third-party libraries:

.. code-block:: cmake

   target_link_libraries(myexe -lsomesystemlib)

If the library does not exist, the end result is a cryptic linker error. See
also `Jussi Pakkanen's blog <http://voices.canonical.com/jussi.pakkanen/2013/03/26/a-list-of-common-cmake-antipatterns/>`_
You will first need to find the library, using the macro
`find_library <http://www.cmake.org/cmake/help/v3.0/command/find_library.html>`_,
and then use the ``target_link_libraries`` command.
