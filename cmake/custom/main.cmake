configure_file (
  "${CMAKE_SOURCE_DIR}/config.h.in"
  "${CMAKE_BINARY_DIR}/config.h"
  )

# FIXME Remove
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set_property(DIRECTORY PROPERTY EP_BASE ${CMAKE_BINARY_DIR}/subprojects)

set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)
message(STATUS "${PROJECT_NAME} staged install: ${STAGED_INSTALL_PREFIX}")

add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(pilot)
