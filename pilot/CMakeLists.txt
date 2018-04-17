link_directories(${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mrchem.cpp)
    add_executable(mrchem-pilot.x mrchem.cpp)
    add_dependencies(mrchem-pilot.x
        mrchem
        submodule-mrcpp
        )
    target_include_directories(mrchem-pilot.x PRIVATE
        ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
        ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/chemistry
        ${CMAKE_SOURCE_DIR}/src/properties
        ${CMAKE_SOURCE_DIR}/src/analyticfunctions
        ${CMAKE_SOURCE_DIR}/src/qmfunctions
        ${CMAKE_SOURCE_DIR}/src/qmoperators
        ${CMAKE_SOURCE_DIR}/src/qmoperators/one_electron
        ${CMAKE_SOURCE_DIR}/src/qmoperators/two_electron
        )
    target_link_libraries(mrchem-pilot.x
        mrchem
        mrcpp
        getkw
        )

    set (MRCHEM_DRIVER_DEBUG ON)
    set (MRCHEM_EXECUTABLE ${CMAKE_BINARY_DIR}/pilot/mrchem-pilot.x)

    configure_file(${CMAKE_SOURCE_DIR}/src/mrchem.in ${CMAKE_BINARY_DIR}/pilot/mrchem-pilot @ONLY)
    execute_process(COMMAND chmod 755 ${CMAKE_BINARY_DIR}/pilot/mrchem-pilot OUTPUT_QUIET)
endif()

if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/mrchem.inp)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mrchem.inp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

