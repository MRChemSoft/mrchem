target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/qmfunctions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QMFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Orbital.cpp
    )

target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
