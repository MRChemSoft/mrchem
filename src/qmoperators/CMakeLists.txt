target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/RankZeroTensorOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QMPotential.cpp
    )

target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )

add_subdirectory(one_electron)
add_subdirectory(two_electron)
