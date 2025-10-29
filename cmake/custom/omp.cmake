#.rst:
#
# Enables OpenMP support.
# This was adapted from Autocmake
#
# Variables used::
#
#   ENABLE_OPENMP
#
# autocmake.yml configuration::
#
#   docopt: "--omp Enable OpenMP parallelization [default: False]."
#   define: "'-DENABLE_OPENMP={0}'.format(arguments['--omp'])"

option(ENABLE_OPENMP "Enable OpenMP parallelization" ON)

if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()
