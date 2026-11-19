option(MRCPP_FIND_BEHAVIOUR "Find MRCPP behavior" "default")
cpm_set_find_behaviour(${MRCPP_FIND_BEHAVIOUR})
set(CMAKE_BUILD_TYPE Release)
set(Eigen3_DIR ${eigen3_BINARY_DIR})
CPMAddPackage(
  NAME MRCPP
  VERSION 1
  GITHUB_REPOSITORY MRChemSoft/mrcpp
  GIT_TAG f73b863d3af98b9bae45d62378a7618285b1dfb1

  FIND_PACKAGE_ARGUMENTS "CONFIG NO_CMAKE_PATH NO_CMAKE_PACKAGE_REGISTRY NO_CMAKE_SYSTEM_PACKAGE_REGISTRY"
  OPTIONS
  "ENABLE_OPENMP ${ENABLE_OPENMP}"
  "ENABLE_MPI ${ENABLE_MPI}"
  "PYTHON_INTERPRETER ${Python_EXECUTABLE}"
  "ENABLE_TESTS OFF"
  "ENABLE_EXAMPLES OFF"
  )

if(NOT mrcpp_ADDED)
  # found locally — MRCPP was pre-built independently of MRChem's own
  # ENABLE_OPENMP/ENABLE_MPI, so check the parallel configs actually match
  get_target_property(MRCPP_HAS_OMP MRCPP::mrcpp MRCPP_HAS_OMP)
  if(ENABLE_OPENMP AND NOT MRCPP_HAS_OMP)
    message(WARNING
      "You are building MRChem with OpenMP, while using a non-OpenMP version of MRCPP!\
         We recommend rebuilding MRCPP with OpenMP support."
      )
  endif()

  get_target_property(MRCPP_HAS_MPI MRCPP::mrcpp MRCPP_HAS_MPI)
  if(ENABLE_MPI AND NOT MRCPP_HAS_MPI)
    message(FATAL_ERROR
      "You cannot build MRChem with MPI and link against a non-MPI version of MRCPP!\
         Rebuild MRCPP with MPI support or disable it for MRChem."
      )
  endif()
endif()
