configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )

set(PYTHON_SITE_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)

# Fetch dependencies
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_getkw.cmake)
file(
  COPY
    ${getkw_PYMOD}/getkw.py
    ${getkw_PYMOD}/pyparsing.py
  DESTINATION
    ${PROJECT_BINARY_DIR}/${PYTHON_SITE_INSTALL_DIR}
  )
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_xcfun.cmake)
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_eigen3.cmake)
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_mrcpp.cmake)

add_subdirectory(src)
add_subdirectory(pilot)
