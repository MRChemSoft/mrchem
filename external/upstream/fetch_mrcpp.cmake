find_package(MRCPP CONFIG QUIET
  NO_CMAKE_PATH
  NO_CMAKE_PACKAGE_REGISTRY
  NO_CMAKE_SYSTEM_PACKAGE_REGISTRY
  )
if(TARGET MRCPP::mrcpp)
  get_property(_loc TARGET MRCPP::mrcpp PROPERTY LOCATION)
  message(STATUS "Found MRCPP: ${_loc} (found version ${MRCPP_VERSION})")

  # check that the parallel configurations of MRChem and MRCPP are compatible
  # these checks are only needed when picking up an installed library:
  # if we build it ourselves, then the parallel configuration for MRCPP will follow that of MRChem
  #
  # 1. OMP MRChem + non-OMP MRCPP is not a great idea, but it's not problematic.
  #    We just emit a warning.
  get_target_property(MRCPP_HAS_OMP MRCPP::mrcpp MRCPP_HAS_OMP)
  if(ENABLE_OPENMP AND NOT MRCPP_HAS_OMP)
    message(WARNING
      "You are building MRChem with OpenMP, while using a non-OpenMP version of MRCPP!\
         We recommend rebuilding MRCPP with OpenMP support."
      )
  endif()

  # 1. MPI MRChem + non-MPI MRCPP will lead to runtime failures.
  #    Fail configuration with a fatal error.
  get_target_property(MRCPP_HAS_MPI MRCPP::mrcpp MRCPP_HAS_MPI)
  if(ENABLE_MPI AND NOT MRCPP_HAS_MPI)
    message(FATAL_ERROR
      "You cannot build MRChem with MPI and link against a non-MPI version of MRCPP!\
         Rebuild MRCPP with MPI support or disable it for MRChem."
      )
  endif()
else()
  message(STATUS "Suitable MRCPP could not be located. Fetching and building!")
  include(FetchContent)

  FetchContent_Declare(mrcpp_sources
    QUIET
    GIT_REPOSITORY
      https://github.com/MRChemSoft/mrcpp.git
    GIT_TAG
      0e5b32003eed8e1abb22a44ba4e1a2ce9fd9ddb1
  )

  FetchContent_GetProperties(mrcpp_sources)

  set(CMAKE_BUILD_TYPE Release)
  set(ENABLE_OPENMP ${ENABLE_OPENMP})
  set(ENABLE_MPI ${ENABLE_MPI})
  set(Eigen3_DIR ${eigen3_sources_BINARY_DIR})
  set(PYTHON_INTERPRETER ${Python_EXECUTABLE})
  set(ENABLE_TESTS ON CACHE BOOL "" FORCE)
  set(ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)

  if(NOT mrcpp_sources_POPULATED)
    FetchContent_Populate(mrcpp_sources)

    add_subdirectory(
      ${mrcpp_sources_SOURCE_DIR}
      ${mrcpp_sources_BINARY_DIR}
      )
  endif()
endif()
