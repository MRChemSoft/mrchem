target_sources(unit-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/identity_operator.cpp
    )

add_Catch_test(identity_operator "identity_operator")
