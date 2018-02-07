include(DefaultPathsMRChem)
if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mrchem.cpp)
    add_executable(mrchem-pilot.x mrchem.cpp)
    add_dependencies(mrchem-pilot.x external-modules)
    target_link_libraries(mrchem-pilot.x mrcpp)
endif()

