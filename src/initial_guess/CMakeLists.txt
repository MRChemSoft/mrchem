link_directories(${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

# Install sources in lib
target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/core.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gto.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sad.cpp
)

#####################################
# Install initial guess executables #
#####################################

# Install core guess executable in bin
add_executable(core-guess core_guess.cpp)
#add_dependencies(core-guess external-modules)
target_include_directories(core-guess PRIVATE
    ${CMAKE_BINARY_DIR}
    ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP
    ${CMAKE_SOURCE_DIR}/src
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
target_link_libraries(core-guess
    mrchem
    mrcpp
    getkw
    )
install(TARGETS core-guess
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install SAD guess executable in bin
add_executable(sad-guess sad_guess.cpp)
#add_dependencies(sad-guess external-modules)
target_include_directories(sad-guess PRIVATE
    ${CMAKE_BINARY_DIR}
    ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP
    ${CMAKE_SOURCE_DIR}/src
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
    ${CMAKE_SOURCE_DIR}/src/initial_guess/sad
    )
target_link_libraries(sad-guess
    mrchem
    mrcpp
    getkw
    )
install(TARGETS sad-guess
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install GTO guess executable in bin
add_executable(gto-guess gto_guess.cpp)
#add_dependencies(gto-guess external-modules)
target_include_directories(gto-guess PRIVATE
    ${CMAKE_BINARY_DIR}
    ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP
    ${CMAKE_SOURCE_DIR}/src
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
target_link_libraries(gto-guess
    mrchem
    mrcpp
    getkw
    )
install(TARGETS gto-guess
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)


# Install driver for LSDalton calculations on stallo
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lsdalton/lsdalton-guess.in lsdalton-guess)
install(
  PROGRAMS
    ${CMAKE_CURRENT_BINARY_DIR}/lsdalton-guess
  DESTINATION
    ${CMAKE_INSTALL_BINDIR}
  )

install(
  FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/lsdalton/lsdalton.py
  DESTINATION
    ${PYTHON_SITE_INSTALL_DIR}
  )

