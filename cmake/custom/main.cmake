include(GNUInstallDirs)

configure_file (
    "${CMAKE_SOURCE_DIR}/config.h.in"
    "${CMAKE_BINARY_DIR}/config.h"
    )

set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(${CMAKE_BINARY_DIR})

add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(pilot)
