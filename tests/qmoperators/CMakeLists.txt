target_sources(mrchem-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/identity_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/position_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/momentum_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/kinetic_operator.cpp
    )

add_catch_test(identity_operator "identity_operator")
add_catch_test(position_operator "position_operator")
add_catch_test(momentum_operator "momentum_operator")
add_catch_test(kinetic_operator "kinetic_operator")
