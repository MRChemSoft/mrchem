find_package(Eigen3 3.3 CONFIG REQUIRED)
if(TARGET Eigen3::Eigen)
  message(STATUS "Using Eigen3: ${EIGEN3_ROOT_DIR} (version ${Eigen3_VERSION})")
endif()
