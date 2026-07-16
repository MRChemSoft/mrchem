if(EIGEN3_FIND_BEHAVIOUR STREQUAL "default" OR EIGEN3_FIND_BEHAVIOUR STREQUAL "onlylocal")
  find_package(Eigen3 3.4 CONFIG QUIET
    NO_CMAKE_PATH
    NO_CMAKE_PACKAGE_REGISTRY
    NO_CMAKE_SYSTEM_PACKAGE_REGISTRY
    )
endif()
if(TARGET Eigen3::Eigen)
  message(STATUS "Found Eigen3: ${EIGEN3_ROOT_DIR} (version ${Eigen3_VERSION})")
elseif(EIGEN3_FIND_BEHAVIOUR STREQUAL "default" OR EIGEN3_FIND_BEHAVIOUR STREQUAL "onlyfetch")
  message(STATUS "Suitable Eigen3 could not be located. Fetching and building!")
  include(FetchContent)

  FetchContent_Declare(eigen3
    QUIET
    URL
      https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )

  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(eigen3)
else()
  message(FATAL_ERROR "No suitable Eigen3 found or fetched. Aborting setup!")
endif()
