target_sources(mrchem-tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/density.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/orbital.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/orbital_vector.cpp
  )

add_Catch_test(
  NAME density
  LABELS "density"
  )

add_Catch_test(
  NAME orbital
  LABELS "orbital"
  )

add_Catch_test(
  NAME orbital_vector
  LABELS "orbital_vector"
  )
