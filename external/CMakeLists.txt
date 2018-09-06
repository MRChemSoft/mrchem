add_custom_target(git-submodule-init
    COMMAND git submodule init
    COMMAND git submodule update
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

find_package(MRCPP CONFIG QUIET)
if(TARGET MRCPP::mrcpp)
  get_property(_loc TARGET MRCPP::mrcpp PROPERTY LOCATION)
  message(STATUS "Found MRCPP: ${_loc} (found version ${MRCPP_VERSION})")
  add_library(mrcpp_external INTERFACE)  # dummy
else()
  message(STATUS "Suitable MRCPP could not be located. Building the local git submodule!")
  include(ExternalProject)
  ExternalProject_Add(mrcpp_external
    DOWNLOAD_COMMAND
      ""
    SOURCE_DIR
      ${CMAKE_CURRENT_SOURCE_DIR}/mrcpp
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DENABLE_OPENMP=${ENABLE_OPENMP}
      -DENABLE_MPI=${ENABLE_MPI}
      -DEigen3_DIR=${Eigen3_DIR}
      -DENABLE_TESTS=TRUE
      -DENABLE_EXAMPLES=TRUE
    CMAKE_CACHE_ARGS
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    TEST_BEFORE_INSTALL
      1
    )
  set(MRCPP_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/MRCPP CACHE PATH "Path to internally built MRCPPConfig.cmake" FORCE)
  # This is very annoying, because it updates the submodule unexpectedly
  #add_dependencies(mrcpp_external git-submodule-init)
endif()

find_package(getkw CONFIG QUIET)
if(TARGET getkw::getkw)
  get_property(_loc TARGET getkw::getkw PROPERTY LOCATION)
  message(STATUS "Found getkw: ${_loc} (found version ${getkw_VERSION})")
  add_library(getkw_external INTERFACE)  # dummy
else()
  message(STATUS "Suitable getkw could not be located. Building the local git submodule!")
  include(ExternalProject)
  ExternalProject_Add(getkw_external
    DOWNLOAD_COMMAND
      ""
    SOURCE_DIR
      ${CMAKE_CURRENT_SOURCE_DIR}/libgetkw
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    CMAKE_CACHE_ARGS
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    TEST_BEFORE_INSTALL
      1
    )
  set(getkw_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/getkw CACHE PATH "Path to internally built getkwConfig.cmake" FORCE)
  # This is very annoying, because it updates the submodule unexpectedly
  #add_dependencies(getkw_external git-submodule-init)
endif()

find_package(XCFun CONFIG QUIET)
if(TARGET XCFun::xcfun)
  get_property(_loc TARGET XCFun::xcfun PROPERTY LOCATION)
  message(STATUS "Found XCFun: ${_loc} (found version ${XCFun_VERSION})")
  add_library(xcfun_external INTERFACE)  # dummy
else()
  # FIXME XCFun needs the C compiler for the moment being
  #       Remove when updating to a latest and greatest that removes this annoyance
  enable_language(C)
  message(STATUS "Suitable XCFun could not be located. Building the local git submodule!")
  include(ExternalProject)
  ExternalProject_Add(xcfun_external
      DOWNLOAD_COMMAND
        ""
      SOURCE_DIR
        ${CMAKE_CURRENT_SOURCE_DIR}/xcfun
      CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DENABLE_FC_SUPPORT=FALSE
        -DENABLE_TESTALL=TRUE
        -DENABLE_OPENMP=${ENABLE_OPENMP}
      CMAKE_CACHE_ARGS
        -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
        -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      TEST_BEFORE_INSTALL
        1
      )
  set(XCFun_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/XCFun CACHE PATH "Path to internally built XCFunConfig.cmake" FORCE)
  # This is very annoying, because it updates the submodule unexpectedly
  #add_dependencies(xcfun_external git-submodule-init)
endif()

install(
  DIRECTORY
    ${STAGED_INSTALL_PREFIX}/
  DESTINATION
    .
  USE_SOURCE_PERMISSIONS
  PATTERN
    "*.pyc" EXCLUDE
  )
