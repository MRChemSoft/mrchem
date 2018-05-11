link_directories(${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

# Install sources in lib
target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/hydrogen_guess.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gto_guess.cpp
)

# Install hydrogen guess executable in bin
add_executable(hydrogen-guess hydrogen.cpp)
add_dependencies(hydrogen-guess external-modules)
target_include_directories(hydrogen-guess PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP
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
target_link_libraries(hydrogen-guess
    mrchem
    mrcpp
    getkw
    )
install(TARGETS hydrogen-guess
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install gto guess executable in bin
add_executable(gto-guess gto.cpp)
add_dependencies(gto-guess external-modules)
target_include_directories(gto-guess PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP
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
