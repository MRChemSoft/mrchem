set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage${CMAKE_INSTALL_PREFIX})

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
      DOWNLOAD_COMMAND ""
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/mrcpp
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                 -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DENABLE_OPENMP=${ENABLE_OPENMP}
                 -DENABLE_MPI=${ENABLE_MPI}
                 -DENABLE_TESTS=FALSE
                 -DENABLE_EXAMPLES=FALSE
                 -DSTATIC_LIBRARY_ONLY=TRUE 
      CMAKE_CACHE_ARGS
                 -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      INSTALL_COMMAND
                 DESTDIR=${CMAKE_BINARY_DIR}/stage ${CMAKE_MAKE_PROGRAM} install
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
      DOWNLOAD_COMMAND ""
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libgetkw
      CMAKE_ARGS  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                  -DCMAKE_BUILD_TYPE=Release
                  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      CMAKE_CACHE_ARGS
                  -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                  -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      INSTALL_COMMAND
                  DESTDIR=${CMAKE_BINARY_DIR}/stage ${CMAKE_MAKE_PROGRAM} install
      )
    # This is only valid when updating to the latest and greatest Getkw, which is WIP
    set(getkw_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/getkw CACHE PATH "Path to internally built getkwConfig.cmake" FORCE)
    # We set these by hand to pass them down to the core project
    # FIXME Remove as soon as we start using the latest Getkw
    set(getkw_INCLUDEDIR ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH "Path to getkw headers" FORCE)
    set(getkw_LIBRARIES ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBRAYDIR}/lib/libgetkw.a CACHE PATH "Path to getkw library" FORCE)
    # This is very annoying, because it updates the submodule unexpectedly
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
      DOWNLOAD_COMMAND ""
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/xcfun
      CMAKE_ARGS  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                  -DCMAKE_BUILD_TYPE=Release
                  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                  -DENABLE_FC_SUPPORT=FALSE
                  -DSTATIC_LIBRARY_ONLY=TRUE
                  -DENABLE_TESTALL=FALSE
                  -DENABLE_OPENMP=${ENABLE_OPENMP}
      CMAKE_CACHE_ARGS
                  -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                  -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      INSTALL_COMMAND
                  DESTDIR=${CMAKE_BINARY_DIR}/stage ${CMAKE_MAKE_PROGRAM} install
      )
  # This is only valid when updating to the latest and greatest XCFun
  # which I am not doing because I don't want conflicts with Luca's PR
  set(XCFun_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/XCFun CACHE PATH "Path to internally built XCFunConfig.cmake" FORCE)
  # We set these by hand to pass them down to the core project
  # FIXME Remove as soon as we start using the latest XCFun
  set(XCFun_INCLUDEDIR ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH "Path to XCFun headers" FORCE)
  set(XCFun_LIBRARIES ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBRAYDIR}/lib/libxcfun.a CACHE PATH "Path to XCFun library" FORCE)
  # This is very annoying, because it updates the submodule unexpectedly
  #add_dependencies(xcfun_external git-submodule-init)
endif()

install(
    DIRECTORY
    ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/
    DESTINATION
    ${CMAKE_INSTALL_PREFIX}
    USE_SOURCE_PERMISSIONS
    )
