link_directories(${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

set (MRCHEM_EXECUTABLE ${CMAKE_BINARY_DIR}/pilot/mrchem-pilot.x)
set (MRCHEM_DRIVER_DEBUG ON)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mrchem.cpp)
    add_executable(mrchem-pilot.x mrchem.cpp)
    add_dependencies(mrchem-pilot.x
        mrchem
        submodule-mrcpp
        )
    target_link_libraries(mrchem-pilot.x
        mrchem
        mrcpp
        getkw
        )
    target_include_directories(mrchem-pilot.x PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/qmfunctions
        ${CMAKE_SOURCE_DIR}/src/qmoperators
        ${CMAKE_SOURCE_DIR}/external/mrcpp/${CMAKE_INSTALL_INCLUDEDIR}
        ${CMAKE_SOURCE_DIR}/external/libgetkw/C++
        )
    configure_file(${CMAKE_SOURCE_DIR}/src/mrchem.in ${CMAKE_BINARY_DIR}/pilot/mrchem)
    execute_process(COMMAND chmod 755 ${CMAKE_BINARY_DIR}/pilot/mrchem OUTPUT_QUIET)
endif()

if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/mrchem.inp)
    file(COPY mrchem.inp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

