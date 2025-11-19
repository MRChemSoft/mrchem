# 1) Prefer a system-provided LibXC with a proper CMake config package
find_package(Libxc QUIET CONFIG)

if(Libxc_FOUND)
  # LibXC should export Libxc::xc; if only 'xc' exists, provide a stable alias
  if(NOT TARGET Libxc::xc AND TARGET xc)
    add_library(Libxc::xc ALIAS xc)
  endif()
  return()  
endif()

include(FetchContent)

# Build settings
set(BUILD_SHARED_LIBS ON  CACHE BOOL "Build LibXC shared libs")
set(BUILD_TESTING     OFF CACHE BOOL "Build LibXC tests")
set(ENABLE_FORTRAN    OFF CACHE BOOL "Build LibXC Fortran bindings")
set(DISABLE_FXC       ON CACHE BOOL "Disable 2nd derivatives (Fxc)")
set(DISABLE_KXC       ON CACHE BOOL "Disable 3rd derivatives (Kxc)")
set(DISABLE_LXC       ON CACHE BOOL "Disable 4th derivatives (Lxc)")  # keep 4th derivs

FetchContent_Declare(libxc_sources
  QUIET
  GIT_REPOSITORY
    https://gitlab.com/libxc/libxc
  GIT_TAG
  7bd5bb41415968db94c499a2f093309c9a2dcf53
)

FetchContent_MakeAvailable(libxc_sources)

# LibXC's CMake exports a target named 'xc'
if(NOT TARGET xc)
  message(FATAL_ERROR "LibXC target 'xc' not found after FetchContent")
endif()

add_library(Libxc::xc ALIAS xc)

get_target_property(_xc_idirs xc INTERFACE_INCLUDE_DIRECTORIES)
if(NOT _xc_idirs)
  # These two dirs contain xc_version.h (build) and xc.h (source)
  target_include_directories(xc PUBLIC
    ${libxc_sources_BINARY_DIR}/include
    ${libxc_sources_SOURCE_DIR}/src
  )
endif()