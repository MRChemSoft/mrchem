
set(HIRSHFELD_SOURCE_DIR ${PROJECT_SOURCE_DIR}/share/hirshfeld CACHE STRING "Path to azora potentials")
set(HIRSHFELD_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/hirshfeld)
install(
  DIRECTORY share/hirshfeld
  DESTINATION share/${PROJECT_NAME}
  )