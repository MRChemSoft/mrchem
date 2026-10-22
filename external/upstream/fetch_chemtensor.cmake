
option(ENABLE_CHEMTENSOR "Enable support for the chemtensor library" OFF)

if(ENABLE_CHEMTENSOR)
  CPMAddPackage(
    NAME Chemtensor
    VERSION 0
    GITHUB_REPOSITORY qc-tum/chemtensor
    GIT_TAG 18c269a5f1c4844c6d68188baa515aa4c9329611
    OPTIONS
    "BUILD_SHARED_LIBS OFF"
    "CMAKE_BUILD_TYPE Debug"
  )
endif()
