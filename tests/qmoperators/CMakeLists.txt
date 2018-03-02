target_sources(mrchem-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/identity_operator.cpp
    )

add_catch_test(identity_operator "identity_operator")
