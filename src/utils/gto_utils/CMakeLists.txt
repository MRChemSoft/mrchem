target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/AOBasis.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AOContraction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Intgrl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OrbitalExp.cpp
)

target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
