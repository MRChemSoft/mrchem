
link_directories(${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})


add_library(mrchem STATIC
    parallel.cpp
    mrenv.cpp
    )
add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)

add_dependencies(mrchem external-modules)
target_link_libraries(mrchem
    mrcpp
    )
target_include_directories(mrchem PRIVATE
    ${CMAKE_SOURCE_DIR}/src/qmfunctions
    ${CMAKE_SOURCE_DIR}/src/qmoperators
    ${CMAKE_SOURCE_DIR}/external/mrcpp/${CMAKE_INSTALL_INCLUDEDIR}
    ${CMAKE_SOURCE_DIR}/external/libgetkw/C++
    )
install(TARGETS mrchem ARCHIVE DESTINATION lib)


add_executable(mrchem.x mrchem.cpp)
add_dependencies(mrchem.x
    mrchem
    )
target_link_libraries(mrchem.x
    mrchem
    mrcpp
    getkw
    )
target_include_directories(mrchem.x PRIVATE
    ${CMAKE_SOURCE_DIR}/src/qmfunctions
    ${CMAKE_SOURCE_DIR}/src/qmoperators
    ${CMAKE_SOURCE_DIR}/external/mrcpp/${CMAKE_INSTALL_INCLUDEDIR}
    ${CMAKE_SOURCE_DIR}/external/libgetkw/C++
    )
set (MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/bin/mrchem.x)
option (MRCHEM_DRIVER_DEBUG "Debug the driver (extmod)" OFF)
install(TARGETS mrchem.x RUNTIME DESTINATION bin)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mrchem.in mrchem)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/mrchem DESTINATION bin)

