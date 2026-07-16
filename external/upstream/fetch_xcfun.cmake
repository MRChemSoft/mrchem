if(XCFUN_FIND_BEHAVIOUR STREQUAL "default" OR XCFUN_FIND_BEHAVIOUR STREQUAL "onlylocal")
  find_package(XCFun 2.1 CONFIG QUIET
    NO_CMAKE_PATH
    NO_CMAKE_PACKAGE_REGISTRY
    NO_CMAKE_SYSTEM_PACKAGE_REGISTRY
    )
endif()
set(XCFUN_LIBRARY XCFun::xcfun)
set(ENABLE_XCFUN ON)
if(TARGET XCFun::xcfun)
  get_property(_loc TARGET XCFun::xcfun PROPERTY LOCATION)
  message(STATUS "Found XCFun: ${_loc} (found version ${XCFun_VERSION})")
elseif(XCFUN_FIND_BEHAVIOUR STREQUAL "default" OR XCFUN_FIND_BEHAVIOUR STREQUAL "onlyfetch")
  message(STATUS "Suitable XCFun could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(xcfun_sources
    QUIET
    URL
      https://github.com/dftlibs/xcfun/archive/v2.1.1.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )

  set(CMAKE_BUILD_TYPE Release)
  set(ENABLE_TESTALL FALSE CACHE BOOL "" FORCE)
  set(XCFUN_MAX_ORDER 3)
  set(XCFUN_PYTHON_INTERFACE FALSE CACHE BOOL "" FORCE)

  if(XCFUN_OLD_PBE)
    message(STATUS "Compiling XCFun with old PBE parameters (different from LibXC)")
  else()
    message(STATUS "Compiling XCFun with new PBE parameters (same as LibXC)")
    add_compile_definitions(XCFUN_REF_PBEX_MU)
  endif()

  FetchContent_MakeAvailable(xcfun_sources)

  if(TARGET xcfun)
    # Suppress all warnings from XCFun's own compilation (-w overrides any -Wall/-Wextra from CMAKE_CXX_FLAGS)
    target_compile_options(xcfun PRIVATE -w)
  endif()
else()
  message(STATUS "Disabling XCFun support!")
  add_compile_definitions(DISABLE_XCFUN)
  set(XCFUN_LIBRARY "")
  set(ENABLE_XCFUN OFF)
endif()
