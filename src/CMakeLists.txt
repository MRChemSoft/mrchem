link_directories(${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

add_library(mrchem STATIC
    parallel.cpp
    mrenv.cpp
    )
add_dependencies(mrchem
    external-modules
    )
add_subdirectory(chemistry)
add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)
target_include_directories(mrchem PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${CMAKE_SOURCE_DIR}/src/chemistry
    ${CMAKE_SOURCE_DIR}/src/properties
    ${CMAKE_SOURCE_DIR}/src/analyticfunctions
    ${CMAKE_SOURCE_DIR}/src/qmfunctions
    ${CMAKE_SOURCE_DIR}/src/qmoperators
    )
target_link_libraries(mrchem
    mrcpp
    )
install(TARGETS mrchem ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/mrchem)


add_executable(mrchem.x
    mrchem.cpp
    )
add_dependencies(mrchem.x
    mrchem
    )
target_include_directories(mrchem.x PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${CMAKE_SOURCE_DIR}/src/chemistry
    ${CMAKE_SOURCE_DIR}/src/properties
    ${CMAKE_SOURCE_DIR}/src/analyticfunctions
    ${CMAKE_SOURCE_DIR}/src/qmfunctions
    ${CMAKE_SOURCE_DIR}/src/qmoperators
    )
target_link_libraries(mrchem.x
    mrchem
    mrcpp
    getkw
    )
install(TARGETS mrchem.x RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})


option (MRCHEM_DRIVER_DEBUG "Debug the driver (extmod)" OFF)
set (MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/mrchem.x)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mrchem.in ${CMAKE_CURRENT_BINARY_DIR}/mrchem)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/mrchem DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

