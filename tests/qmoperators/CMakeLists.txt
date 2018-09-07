target_sources(mrchem-tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/identity_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/position_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/momentum_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/kinetic_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nuclear_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/electric_field_operator.cpp
  )

add_Catch_test(
  NAME identity_operator
  LABELS "identity_operator"
  )

add_Catch_test(
  NAME position_operator
  LABELS "position_operator"
  )

add_Catch_test(
  NAME momentum_operator
  LABELS "momentum_operator"
  )

add_Catch_test(
  NAME kinetic_operator
  LABELS "kinetic_operator"
  )

add_Catch_test(
  NAME nuclear_operator
  LABELS "nuclear_operator"
  )

add_Catch_test(
  NAME electric_field_operator
  LABELS "electric_field_operator"
  )
