link_directories(${STAGED_INSTALL_PREFIX}/lib)
link_directories(${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

add_library(mrchem STATIC
    SCFDriver.cpp
    parallel.cpp
    mrenv.cpp
    )
  #add_dependencies(mrchem
  #    external-modules
  #    )
add_subdirectory(chemistry)
add_subdirectory(analyticfunctions)
add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)
add_subdirectory(mrdft)
add_subdirectory(scf_solver)
add_subdirectory(initial_guess)
add_subdirectory(utils)
target_include_directories(mrchem PRIVATE
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
install(TARGETS mrchem.x RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})


option (MRCHEM_DRIVER_DEBUG "Debug the driver (extmod)" OFF)
set (MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/mrchem.x)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mrchem.in ${CMAKE_CURRENT_BINARY_DIR}/mrchem)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/mrchem DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
