set(AZORA_POTENTIALS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/share/azora_potentials CACHE STRING "Path to azora potentials")
set(AZORA_POTENTIALS_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/azora_potentials)
install(
  DIRECTORY share/azora_potentials
  DESTINATION share/${PROJECT_NAME}
  )

