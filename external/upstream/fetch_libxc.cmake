if(LIBXC_FIND_BEHAVIOUR STREQUAL "default" OR LIBXC_FIND_BEHAVIOUR STREQUAL "onlylocal")
  find_package(Libxc 7.0 CONFIG QUIET
    NO_CMAKE_PATH
    NO_CMAKE_PACKAGE_REGISTRY
    NO_CMAKE_SYSTEM_PACKAGE_REGISTRY
    )
endif()
if(TARGET Libxc::xc)
  get_target_property(_loc Libxc::xc LOCATION)
  message(STATUS "Found Libxc: ${_loc} (found version ${Libxc_VERSION})")
elseif(LIBXC_FIND_BEHAVIOUR STREQUAL "default" OR LIBXC_FIND_BEHAVIOUR STREQUAL "onlyfetch")
  message(STATUS "Suitable LibXC could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(libxc_sources
    QUIET
    URL
      https://gitlab.com/libxc/libxc/-/archive/7.0.0/libxc-7.0.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )

  set(BUILD_TESTING     OFF CACHE BOOL "Build LibXC tests"          FORCE)
  set(ENABLE_TESTS      OFF CACHE BOOL "Enable LibXC tests"         FORCE)
  set(BUILD_SHARED_LIBS ON  CACHE BOOL "Build LibXC shared libs"    FORCE)
  set(ENABLE_FORTRAN    OFF CACHE BOOL "Build LibXC Fortran bindings" FORCE)
  set(DISABLE_FXC       ON  CACHE BOOL "Disable 2nd derivatives (Fxc)" FORCE)
  set(DISABLE_KXC       ON  CACHE BOOL "Disable 3rd derivatives (Kxc)" FORCE)
  set(DISABLE_LXC       ON  CACHE BOOL "Disable 4th derivatives (Lxc)" FORCE)

  # LibXC installs headers flat (no own subdirectory), so we redirect its
  # install include dir to a namespaced folder and restore afterward.
  set(_saved_install_includedir "${CMAKE_INSTALL_INCLUDEDIR}")
  set(CMAKE_INSTALL_INCLUDEDIR "include/LibXC" CACHE STRING "" FORCE)
  FetchContent_MakeAvailable(libxc_sources)
  set(CMAKE_INSTALL_INCLUDEDIR "${_saved_install_includedir}" CACHE STRING "" FORCE)

  if(TARGET xc)
    if(NOT TARGET Libxc::xc)
      add_library(Libxc::xc ALIAS xc)
    endif()

    target_include_directories(xc
      INTERFACE
        $<BUILD_INTERFACE:${libxc_sources_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${libxc_sources_BINARY_DIR}>
        $<INSTALL_INTERFACE:include/LibXC>
      )
  endif()
else()
  message(FATAL_ERROR "No suitable LibXC found or fetched. Aborting setup!")
endif()
