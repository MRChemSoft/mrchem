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
    ${CMAKE_CURRENT_SOURCE_DIR}/operator_composition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/two_el_int.cpp
  )

add_Catch_test(
  NAME identity_operator
  LABELS "identity;operator;short"
  )

add_Catch_test(
  NAME position_operator
  LABELS "position;operator;short"
  )

add_Catch_test(
  NAME momentum_operator
  LABELS "momentum;operator;short"
  )

add_Catch_test(
  NAME kinetic_operator
  LABELS "kinetic;operator;short"
  )

add_Catch_test(
  NAME nuclear_operator
  LABELS "nuclear;operator;short"
  )

add_Catch_test(
  NAME coulomb_operator
  LABELS "coulomb;operator;short"
  )

add_Catch_test(
  NAME coulomb_hessian
  LABELS "coulomb;hessian;short"
  )

add_Catch_test(
  NAME exchange_operator
  LABELS "exchange;operator;short"
  )

add_Catch_test(
  NAME exchange_hessian
  LABELS "exchange;hessian;short"
  )

add_Catch_test(
  NAME xc_operator_lda
  LABELS "xc;operator;lda;short"
  )

add_Catch_test(
  NAME xc_operator_blyp
  LABELS "xc;operator;blyp;short"
  )

add_Catch_test(
  NAME xc_hessian_lda
  LABELS "xc;hessian;lda;short"
  )

add_Catch_test(
  NAME xc_hessian_pbe
  LABELS "xc;hessian;pbe;short"
  )

add_Catch_test(
  NAME electric_field_operator
  LABELS "electric_field;operator;short"
  )

add_Catch_test(
  NAME electric_field_energy
  LABELS "electric_field;energy;short"
  )

add_Catch_test(
  NAME operator_composition
  LABELS "operator;composition;short"
  )

add_Catch_test(
  NAME two_el_int
  LABELS "two-electron;integration;short"
  )
