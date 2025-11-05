find_package(XCFun CONFIG QUIET
  NO_CMAKE_PATH
  NO_CMAKE_PACKAGE_REGISTRY
  NO_CMAKE_SYSTEM_PACKAGE_REGISTRY
  )
if(TARGET XCFun::xcfun)
  get_property(_loc TARGET XCFun::xcfun PROPERTY LOCATION)
  message(STATUS "Found XCFun: ${_loc} (found version ${XCFun_VERSION})")
else()
  message(STATUS "Suitable XCFun could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(xcfun_sources
    QUIET
    URL
      https://github.com/dftlibs/xcfun/archive/v2.1.0.tar.gz
    )

  FetchContent_GetProperties(xcfun_sources)

  set(CMAKE_BUILD_TYPE Release)
  set(ENABLE_TESTALL FALSE CACHE BOOL "")
  set(XCFUN_MAX_ORDER 3)  # TODO Maybe as a user-facing option?
  set(XCFUN_PYTHON_INTERFACE FALSE CACHE BOOL "")

  if(NOT xcfun_sources_POPULATED)
    FetchContent_Populate(xcfun_sources)

    add_subdirectory(
      ${xcfun_sources_SOURCE_DIR}
      ${xcfun_sources_BINARY_DIR}
      )
      message(xcfun_sources_SOURCE_DIR="${xcfun_sources_SOURCE_DIR}")
      message(xcfun_sources_BINARY_DIR="${xcfun_sources_BINARY_DIR}")
  endif()
endif()
