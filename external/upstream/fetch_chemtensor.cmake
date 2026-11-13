
option(ENABLE_CHEMTENSOR "Enable support for the chemtensor library" OFF)

if(ENABLE_CHEMTENSOR)
  CPMAddPackage(
    NAME Chemtensor
    VERSION 0
    GITHUB_REPOSITORY qc-tum/chemtensor
    GIT_TAG b46e1bb75eac11dfd0e14e2fd66f2f058f6f67a0
    OPTIONS
    "BUILD_SHARED_LIBS OFF"
    "ENABLE_CHEMTENSOR_UNIT_TESTS OFF"
    "ENABLE_CHEMTENSOR_PYTHON_MODULE OFF"
    "CMAKE_BUILD_TYPE Release"
  )
endif()
