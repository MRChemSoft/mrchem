include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/external/catch)
include_directories(${CMAKE_BINARY_DIR})

link_directories(${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

add_executable(unit-tests unit_tests.cpp)
add_dependencies(unit-tests submodule-mrcpp)

add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)

target_link_libraries(unit-tests
    mrchem
    mrcpp
    )
    
target_include_directories(unit-tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/qmfunctions
    ${CMAKE_SOURCE_DIR}/src/qmoperators
    ${CMAKE_SOURCE_DIR}/external/mrcpp/${CMAKE_INSTALL_INCLUDEDIR}
    )
