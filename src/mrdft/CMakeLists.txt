target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Factory.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MRDFT.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Functional.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/xc_utils.cpp
    )

add_subdirectory(xc_libraries)
