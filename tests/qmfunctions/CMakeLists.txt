target_sources(mrchem-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/density.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/orbital.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/orbital_vector.cpp
    )

add_catch_test(density "density")
add_catch_test(orbital "orbital")
add_catch_test(orbital_vector "orbital_vector")
