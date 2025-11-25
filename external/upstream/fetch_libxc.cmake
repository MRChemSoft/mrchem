cmake_policy(SET CMP0144 NEW)
find_package(Libxc QUIET CONFIG)

if(TARGET Libxc::xc)
  get_property(_loc TARGET Libxc::xc PROPERTY LOCATION)
  message(STATUS "Libxc::xc: ${_loc} (found version ${Libxc_VERSION})")
else()
  message(STATUS "Suitable LibXC could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(libxc_sources
    QUIET
    GIT_REPOSITORY
      https://gitlab.com/libxc/libxc
    GIT_TAG
    7bd5bb41415968db94c499a2f093309c9a2dcf53
  )

  FetchContent_MakeAvailable(libxc_sources)
  FetchContent_GetProperties(libxc_sources)

  set(BUILD_SHARED_LIBS ON  CACHE BOOL "Build LibXC shared libs")
  set(BUILD_TESTING     OFF CACHE BOOL "Build LibXC tests")
  set(ENABLE_FORTRAN    OFF CACHE BOOL "Build LibXC Fortran bindings")
  set(DISABLE_FXC       ON  CACHE BOOL "Disable 2nd derivatives (Fxc)")
  set(DISABLE_KXC       ON  CACHE BOOL "Disable 3rd derivatives (Kxc)")
  set(DISABLE_LXC       ON  CACHE BOOL "Disable 4th derivatives (Lxc)")

  add_library(Libxc::xc ALIAS xc)

  if(NOT libxc_sources_POPULATED)
    FetchContent_Populate(libxc_sources)
    add_subdirectory(
      ${libxc_sources_SOURCE_DIR}
      ${libxc_sources_BINARY_DIR}
      )
  endif()
endif()
