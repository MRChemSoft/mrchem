add_subdirectory(gto_utils)

target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/mathutils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NonlinearMaximizer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RRMaximizer.cpp
    )
