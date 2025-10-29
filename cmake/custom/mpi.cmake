#.rst:
#
# Enables MPI support.
# This was adapted from Autocmake
#
# Variables used::
#
#   ENABLE_MPI
#
# autocmake.yml configuration::
#
#   docopt: "--mpi Enable MPI parallelization [default: False]."
#   define: "'-DENABLE_MPI={0}'.format(arguments['--mpi'])"

option(ENABLE_MPI "Enable MPI parallelization" ON)

if(ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()
