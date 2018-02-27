include(DefaultPathsMRChem)
set_property(GLOBAL PROPERTY TestSources)
add_subdirectory(mrcpp)

get_property(test_sources GLOBAL PROPERTY TestSources)

include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR}/external/catch)
include_directories(${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories(${CMAKE_SOURCE_DIR}/src/mrcpp/mwanalytic)
include_directories(${CMAKE_SOURCE_DIR}/src/mrcpp/mwcore)
include_directories(${CMAKE_SOURCE_DIR}/src/mrcpp/mwtrees)
include_directories(${CMAKE_SOURCE_DIR}/src/mrcpp/mwoperators)
include_directories(${CMAKE_SOURCE_DIR}/src/mrcpp/mwbuilders)

add_executable(unit_tests.x unit_tests.cpp ${test_sources})

target_link_libraries(unit_tests.x 
    mrcpp 
    mwanalytic
    mwcore
    mwtrees
    mwoperators
    mwbuilders
)
