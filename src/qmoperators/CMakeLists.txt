target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/QMDerivative.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QMIdentity.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QMPotential.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QMSpin.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/qmoperator_utils.cpp
    )

add_subdirectory(one_electron)
add_subdirectory(two_electron)
