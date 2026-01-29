target_sources(mrchem-tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/qmfunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/density.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/orbital.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/orbital_vector.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/spherical_harmonics.cpp
  )

add_Catch_test(
  NAME spherical_harmonics
  LABELS "spherical_harmonics"
  )

add_Catch_test(
  NAME qmfunction
  LABELS "qmfunction"
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
