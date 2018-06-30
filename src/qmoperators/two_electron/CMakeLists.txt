target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/CoulombPotential.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExchangePotential.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FockOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/XCPotential.cpp
    )

target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
