# Use a version range so that newer Eigen releases are accepted: with modern
# Eigen3ConfigVersion.cmake a bare "3.4" means "3.4.x only" (3.4 <= v < 3.5),
# which rejects e.g. Eigen 5.x. The range matches 3.4.0 up to (but not incl.) 6.0.
find_package(Eigen3 3.4...6 CONFIG QUIET
  NO_CMAKE_PATH
  NO_CMAKE_PACKAGE_REGISTRY
  )

if(TARGET Eigen3::Eigen)
  message(STATUS "Using Eigen3: ${EIGEN3_ROOT_DIR} (version ${Eigen3_VERSION})")
else()
  message(STATUS "Suitable Eigen3 could not be located. Fetching and building!")
  include(FetchContent)

  FetchContent_Declare(eigen3
    QUIET
    URL
      https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    )

  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(eigen3)
endif()
