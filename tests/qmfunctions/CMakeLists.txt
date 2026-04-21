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
  LABELS "spherical_harmonics;short"
  )

add_Catch_test(
  NAME qmfunction
  LABELS "qmfunction;short"
  )

add_Catch_test(
  NAME density
  LABELS "density;short"
  )

add_Catch_test(
  NAME orbital
  LABELS "orbital;short"
  )

add_Catch_test(
  NAME orbital_vector
  LABELS "orbital_vector;short"
  )
