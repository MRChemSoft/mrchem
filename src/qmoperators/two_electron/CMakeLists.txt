target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/CoulombPotential.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CoulombPotentialD1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CoulombPotentialD2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExchangePotential.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExchangePotentialD1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExchangePotentialD2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FockOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/XCPotential.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/XCPotentialD1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/XCPotentialD2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ReactionPotential.cpp
    )
