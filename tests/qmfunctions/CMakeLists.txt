target_sources(unit-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/orbital.cpp
    )

add_Catch_test(orbital "orbital")
