target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/chemistry_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PeriodicTable.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Molecule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Cavity.cpp
    )

# copy periodictable.py in mrchem folder (in the build folder)
file(
  COPY
    ${CMAKE_CURRENT_LIST_DIR}/periodictable.py
  DESTINATION
    ${CMAKE_BINARY_DIR}/${PYTHON_SITE_INSTALL_DIR}/mrchem
  )
