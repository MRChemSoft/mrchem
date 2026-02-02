find_package(Libxc QUIET CONFIG)

if(TARGET Libxc::xc)
  get_target_property(_loc Libxc::xc LOCATION)
  message(STATUS "Libxc::xc: ${_loc} (found version ${Libxc_VERSION})")
else()
  message(STATUS "Suitable LibXC could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(libxc_sources
    GIT_REPOSITORY https://gitlab.com/libxc/libxc
    GIT_TAG 7.0.0
  )
  set(CMAKE_INSTALL_INCLUDEDIR "include/" CACHE STRING "" FORCE) # Creates Libxc subdir in install
  set(BUILD_TESTING     OFF CACHE BOOL "Build LibXC tests" FORCE)
  set(ENABLE_TESTS      OFF CACHE BOOL "Enable LibXC tests" FORCE)
  set(BUILD_SHARED_LIBS ON  CACHE BOOL "Build LibXC shared libs" FORCE)
  set(ENABLE_FORTRAN    OFF CACHE BOOL "Build LibXC Fortran bindings" FORCE)
  set(DISABLE_FXC       ON  CACHE BOOL "Disable 2nd derivatives (Fxc)" FORCE)
  set(DISABLE_KXC       ON  CACHE BOOL "Disable 3rd derivatives (Kxc)" FORCE)
  set(DISABLE_LXC       ON  CACHE BOOL "Disable 4th derivatives (Lxc)" FORCE)

  FetchContent_MakeAvailable(libxc_sources)

  if(TARGET xc)
    if(NOT TARGET Libxc::xc)
      add_library(Libxc::xc ALIAS xc)
    endif()

    target_include_directories(xc
      INTERFACE
        $<BUILD_INTERFACE:${libxc_sources_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${libxc_sources_BINARY_DIR}>
        $<INSTALL_INTERFACE:include/Libxc> 
    )
  endif()
endif()