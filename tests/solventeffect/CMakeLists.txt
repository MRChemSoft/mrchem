target_sources(mrchem-tests
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/reaction_operator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cavity_function.cpp
  )

add_Catch_test(
  NAME cavity_function
  LABELS "cavity_function"
  )

add_Catch_test(
  NAME reaction_operator
  LABELS "reaction_operator"
  )
