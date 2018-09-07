link_directories(${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

# RPATH fixing
file(RELATIVE_PATH _rel ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR} ${CMAKE_INSTALL_PREFIX})
if(APPLE)
  set(_rpath "@loader_path/${_rel}")
else()
  set(_rpath "\$ORIGIN/${_rel}")
endif()
file(TO_NATIVE_PATH "${_rpath}/${CMAKE_INSTALL_LIBDIR}" MRChem_RPATH)

add_library(mrchem STATIC
    SCFDriver.cpp
    parallel.cpp
    mrenv.cpp
    )
  #add_dependencies(mrchem
  #    external-modules
  #    )

target_include_directories(mrchem PUBLIC
    ${CMAKE_BINARY_DIR}
    ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP
    ${CMAKE_SOURCE_DIR}/src/analyticfunctions
    ${CMAKE_SOURCE_DIR}/src/chemistry
    ${CMAKE_SOURCE_DIR}/src/initial_guess
    ${CMAKE_SOURCE_DIR}/src/properties
    ${CMAKE_SOURCE_DIR}/src/qmfunctions
    ${CMAKE_SOURCE_DIR}/src/qmoperators
    ${CMAKE_SOURCE_DIR}/src/qmoperators/one_electron
    ${CMAKE_SOURCE_DIR}/src/qmoperators/two_electron
    ${CMAKE_SOURCE_DIR}/src/scf_solver
    ${CMAKE_SOURCE_DIR}/src/utils
    )

add_subdirectory(chemistry)
add_subdirectory(analyticfunctions)
add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)
add_subdirectory(mrdft)
add_subdirectory(scf_solver)
add_subdirectory(initial_guess)
add_subdirectory(utils)
target_link_libraries(mrchem
    PUBLIC
      Eigen3::Eigen
      mrcpp
      xcfun
    )
install(TARGETS mrchem ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/mrchem)


add_executable(mrchem.x
    mrchem.cpp
    )
  #add_dependencies(mrchem.x
  #    mrchem
  #    )
target_include_directories(mrchem.x PRIVATE
    ${CMAKE_BINARY_DIR}
    ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP
    ${CMAKE_SOURCE_DIR}/src/analyticfunctions
    ${CMAKE_SOURCE_DIR}/src/chemistry
    ${CMAKE_SOURCE_DIR}/src/initial_guess
    ${CMAKE_SOURCE_DIR}/src/properties
    ${CMAKE_SOURCE_DIR}/src/qmfunctions
    ${CMAKE_SOURCE_DIR}/src/qmoperators
    ${CMAKE_SOURCE_DIR}/src/qmoperators/one_electron
    ${CMAKE_SOURCE_DIR}/src/qmoperators/two_electron
    ${CMAKE_SOURCE_DIR}/src/scf_solver
    ${CMAKE_SOURCE_DIR}/src/utils
    )
target_link_libraries(mrchem.x
    mrchem
    mrcpp
    getkw
    xcfun
    )

set_target_properties(mrchem.x
  PROPERTIES
    MACOSX_RPATH ON
    SKIP_BUILD_RPATH OFF
    BUILD_WITH_INSTALL_RPATH OFF
    INSTALL_RPATH "${MRChem_RPATH}"
    INSTALL_RPATH_USE_LINK_PATH ON
  )

install(TARGETS mrchem.x RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})


option (MRCHEM_DRIVER_DEBUG "Debug the driver (extmod)" OFF)
set(MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/mrchem.x)
# configure into CMakeFiles directory, because it must be in a different directory
# since file(COPY) does rename files
configure_file(mrchem.in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mrchem @ONLY)
# now copy the temporary into the final destination, setting the permissions
file(
  COPY ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mrchem
  DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
  FILE_PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
  )
install(
  PROGRAMS
    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/mrchem
  DESTINATION
    ${CMAKE_INSTALL_BINDIR}
  )
