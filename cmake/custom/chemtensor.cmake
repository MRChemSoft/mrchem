#.rst:
#
# Enables ChemTensor support.
# This was adapted from Autocmake
#
# Variables used::
#
#   ENABLE_CHEMTENSOR
#
# autocmake.yml configuration::
#
#   docopt: "--chemtensor Enable ChemTensor support [default: False]."
#   define: "'-DENABLE_CHEMTENSOR={0}'.format(arguments['--chemtensor'])"

option(ENABLE_CHEMTENSOR "Enable support for the chemtensor library" OFF)

if(ENABLE_CHEMTENSOR)
  include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_chemtensor.cmake)
endif()
