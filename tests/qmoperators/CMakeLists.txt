target_sources(mrchem-tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/identity_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/position_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/momentum_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/kinetic_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nuclear_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/coulomb_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/coulomb_hessian.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/exchange_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/exchange_hessian.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/xc_operator_lda.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/xc_operator_blyp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/xc_hessian_lda.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/xc_hessian_pbe.cpp
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
  NAME coulomb_operator
  LABELS "coulomb_operator"
  )

add_Catch_test(
  NAME coulomb_hessian
  LABELS "coulomb_hessian"
  )

add_Catch_test(
  NAME exchange_operator
  LABELS "exchange_operator"
  )

add_Catch_test(
  NAME exchange_hessian
  LABELS "exchange_hessian"
  )

add_Catch_test(
  NAME xc_operator_lda
  LABELS "xc_operator_lda"
  )

add_Catch_test(
  NAME xc_operator_blyp
  LABELS "xc_operator_blyp"
  )

add_Catch_test(
  NAME xc_hessian_lda
  LABELS "xc_hessian_lda"
  )

add_Catch_test(
  NAME xc_hessian_pbe
  LABELS "xc_hessian_pbe"
  )

add_Catch_test(
  NAME electric_field_operator
  LABELS "electric_field_operator"
  )

add_Catch_test(
  NAME electric_field_energy
  LABELS "electric_field_energy"
  )
