include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/external/catch)

link_directories(${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

add_executable(mrchem-tests unit_tests.cpp)
#add_dependencies(mrchem-tests submodule-mrcpp)

add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)

target_include_directories(mrchem-tests PRIVATE
    ${CMAKE_BINARY_DIR}
    ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/chemistry
    ${CMAKE_SOURCE_DIR}/src/properties
    ${CMAKE_SOURCE_DIR}/src/analyticfunctions
    ${CMAKE_SOURCE_DIR}/src/qmfunctions
    ${CMAKE_SOURCE_DIR}/src/qmoperators
    ${CMAKE_SOURCE_DIR}/src/qmoperators/one_electron
    ${CMAKE_SOURCE_DIR}/src/qmoperators/two_electron
    )

target_link_libraries(mrchem-tests
    mrchem
    mrcpp
    )

