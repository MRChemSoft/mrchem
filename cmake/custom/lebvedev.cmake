set(LEBVEDEV_SOURCE_DIR ${PROJECT_SOURCE_DIR}/share/lebvedev CACHE STRING "Path to SAD basis and density files")
set(LEBVEDEV_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/lebvedev)
install(
  DIRECTORY share/lebvedev
  DESTINATION share/${PROJECT_NAME}
  )

