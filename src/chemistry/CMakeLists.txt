target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/chem_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PeriodicTable.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Molecule.cpp
    )

target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )

install(
  FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/periodictable.py
  DESTINATION
    ${PYTHON_SITE_INSTALL_DIR}
  )
