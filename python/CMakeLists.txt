find_package(Python3 COMPONENTS Interpreter REQUIRED)

if("${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}" VERSION_LESS 3.6)
    message(FATAL_ERROR "You will need Python >= 3.6 to run MRChem")
endif()

if(NOT DEFINED PYMOD_INSTALL_LIBDIR)
  message(STATUS "Setting (unspecified) option PYMOD_INSTALL_LIBDIR: python")
  set(PYMOD_INSTALL_LIBDIR "python" CACHE STRING "Location under 'lib' to which Python modules are installed" FORCE)
else()
  message(STATUS "Setting option PYMOD_INSTALL_LIBDIR: ${PYMOD_INSTALL_LIBDIR}")
  set(PYMOD_INSTALL_LIBDIR "${PYMOD_INSTALL_LIBDIR}" CACHE STRING "Location under 'lib' to which Python modules are installed" FORCE)
endif()
file(TO_NATIVE_PATH "lib/${PYMOD_INSTALL_LIBDIR}/mrchem" PYMOD_INSTALL_FULLDIR)

# define variable used in mrchem.in
set(MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/mrchem.x)

# configure the mrchem input parser script. This goes into the CMakeFiles directory,
# because it must be in a different directory since file(COPY) only rename files
configure_file(mrchem.in ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mrchem @ONLY)
# now copy the temporary into the final destination, setting the permissions
file(
  COPY
    ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mrchem
  DESTINATION
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
  FILE_PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
  )

# copy mrchem folder under build tree
add_custom_target(mrchem-python
  ALL
    ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_LIST_DIR}/mrchem/*.py ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
  COMMAND
    ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/mrchem/input_parser ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/input_parser
  COMMENT
    "Copy ${CMAKE_CURRENT_LIST_DIR}/mrchem to ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}"
  )
add_dependencies(mrchem.x mrchem-python)

# configure config.py.in
configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/mrchem/config.py.in
  ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/config.py
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
# NOTE the trailing "/" is important!
install(
  DIRECTORY
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}
  FILES_MATCHING
    PATTERN "*.py"
  PATTERN "__pycache__" EXCLUDE
  PATTERN "*.py.in" EXCLUDE
  )
