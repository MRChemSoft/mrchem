find_package(MRCPP CONFIG QUIET
  NO_CMAKE_PATH
  NO_CMAKE_PACKAGE_REGISTRY
  NO_CMAKE_SYSTEM_PACKAGE_REGISTRY
  )
if(TARGET MRCPP::mrcpp)
  get_property(_loc TARGET MRCPP::mrcpp PROPERTY LOCATION)
  message(STATUS "Found MRCPP: ${_loc} (found version ${MRCPP_VERSION})")
else()
  message(STATUS "Suitable MRCPP could not be located. Fetching and building!")
  include(FetchContent)

  FetchContent_Declare(mrcpp_sources
    QUIET
    URL
      https://github.com/MRChemSoft/mrcpp/archive/v1.3.3.tar.gz
    )

  FetchContent_GetProperties(mrcpp_sources)

  set(CMAKE_BUILD_TYPE Release)
  set(ENABLE_OPENMP ${ENABLE_OPENMP})
  set(ENABLE_MPI ${ENABLE_MPI})
  set(Eigen3_DIR ${eigen3_sources_BINARY_DIR})
  set(PYTHON_INTERPRETER ${PYTHON_EXECUTABLE})
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
