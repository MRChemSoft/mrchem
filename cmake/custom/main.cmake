include(GNUInstallDirs)

configure_file (
    "${CMAKE_SOURCE_DIR}/config.h.in"
    "${CMAKE_BINARY_DIR}/config.h"
    )

add_subdirectory(external)

include(ExternalProject)
ExternalProject_Add(${PROJECT_NAME}_core
  DEPENDS
    mrcpp_external
    getkw_external
    xcfun_external
  SOURCE_DIR
    ${PROJECT_SOURCE_DIR}/src
  CMAKE_ARGS
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=11
    -DCMAKE_CXX_EXTENSIONS=OFF
    -DCMAKE_CXX_STANDARD_REQUIRED=ON
    # FIXME Use Eigen3 imported target, requires Eigen 3.3
    #-DEigen3_DIR=${Eigen3_DIR}
    -DEIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR}
    -DMRCPP_DIR=${MRCPP_DIR}
    -DXCFun_DIR=${XCFun_DIR}
    -DGetkw_DIR=${Getkw_DIR}
    # FIXME
    -DGetkw_INCLUDEDIR=${Getkw_INCLUDEDIR}
    # FIXME
    -DGetkw_LIBRARIES=${Getkw_LIBRARIES}
    -DPYTHON_SITE_INSTALL_DIR=${PYTHON_SITE_INSTALL_DIR}
  CMAKE_CACHE_ARGS
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
  BUILD_ALWAYS
    1
  INSTALL_COMMAND
    ""
  )

add_subdirectory(pilot)
