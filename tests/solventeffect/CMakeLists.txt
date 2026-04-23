target_sources(mrchem-tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/reaction_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cavity_function.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PB_solver.cpp
  )

add_Catch_test(
  NAME cavity_function
  LABELS "cavity;function;short"
  )

add_Catch_test(
  NAME reaction_operator
  LABELS "reaction;operator;short"
  )

 add_Catch_test(
   NAME PB_solver
   LABELS "PB;solver;medium"
   )
