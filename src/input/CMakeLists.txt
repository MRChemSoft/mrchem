# define variable used in mrchem.in
set(MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/mrchem.x)

# configure the mrchem input parser script. This goes into the CMakeFiles directory,
# because it must be in a different directory since file(COPY) only rename files
configure_file(mrchem.in ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mrchem @ONLY)

# now copy the temporary into the final destination, setting the permissions
file(
  COPY ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mrchem
  DESTINATION ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
  FILE_PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
  )

# copy mrchem folder to build folder
file(
  COPY
    ${CMAKE_CURRENT_LIST_DIR}/mrchem/
  DESTINATION
    ${CMAKE_BINARY_DIR}/${PYTHON_SITE_INSTALL_DIR}/mrchem
  PATTERN "config.py.in" EXCLUDE
  )
# configure config.py.in
configure_file(
  mrchem/config.py.in
  ${CMAKE_BINARY_DIR}/${PYTHON_SITE_INSTALL_DIR}/mrchem/config.py
  @ONLY
  )

# install the frontend script alongside the mrchem.x executable
install(
  PROGRAMS
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/mrchem
  DESTINATION
    ${CMAKE_INSTALL_BINDIR}
  )

# install the mrchem folder, which includes the (generated) input_parser folder
install(
  DIRECTORY
    ${CMAKE_BINARY_DIR}/${PYTHON_SITE_INSTALL_DIR}/mrchem
  DESTINATION
    ${PYTHON_SITE_INSTALL_DIR}
  FILES_MATCHING
    PATTERN "*.py"
  PATTERN "__pycache__" EXCLUDE
  PATTERN "*.py.in" EXCLUDE
  )
