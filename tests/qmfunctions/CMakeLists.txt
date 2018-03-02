target_sources(mrchem-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/orbital.cpp
    )

add_catch_test(orbital "orbital")
