find_package(MRCPP CONFIG QUIET)
if(TARGET MRCPP::mrcpp)
  get_property(_loc TARGET MRCPP::mrcpp PROPERTY LOCATION)
  message(STATUS "Found MRCPP: ${_loc} (found version ${MRCPP_VERSION})")
else()
  message(STATUS "Suitable MRCPP could not be located. Fetching and building!")
  include(FetchContent)

  FetchContent_Populate(mrcpp_sources
    QUIET
    GIT_REPOSITORY
      https://github.com/MRChemSoft/mrcpp.git
    GIT_TAG
      v1.0.2
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DENABLE_OPENMP=${ENABLE_OPENMP}
      -DENABLE_MPI=${ENABLE_MPI}
      -DEigen3_DIR=${eigen3_sources_BINARY_DIR}
      -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE} # Seems to be ignored...
      -DENABLE_TESTS=TRUE
      -DENABLE_EXAMPLES=FALSE
    CMAKE_CACHE_ARGS
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    )

  add_subdirectory(
    ${mrcpp_sources_SOURCE_DIR}
    ${mrcpp_sources_BINARY_DIR}
    )
endif()
