add_subdirectory(one_electron)
add_subdirectory(two_electron)

target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/RankZeroTensorOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QMOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QMPotential.cpp
    )
