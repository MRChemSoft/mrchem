target_sources(mrchem PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/HydrogenFunction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/NuclearFunction.cpp
  )

target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
