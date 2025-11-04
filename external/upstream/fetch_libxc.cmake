# Always fetch & build LibXC from source
include(FetchContent)

set(LIBXC_VERSION "7.0.0" CACHE STRING "LibXC version to fetch")

# Build settings
set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
set(BUILD_TESTING     OFF CACHE BOOL "" FORCE)
set(ENABLE_FORTRAN    OFF CACHE BOOL "" FORCE)
set(DISABLE_KXC       OFF CACHE BOOL "" FORCE)  # keep 3rd derivs
set(DISABLE_LXC       OFF CACHE BOOL "" FORCE)  # keep 4th derivs

FetchContent_Declare(libxc_sources
  URL "https://gitlab.com/libxc/libxc/-/archive/${LIBXC_VERSION}/libxc-${LIBXC_VERSION}.tar.gz"
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  QUIET
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
