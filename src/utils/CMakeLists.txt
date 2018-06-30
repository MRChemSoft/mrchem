add_subdirectory(gto_utils)

target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/math_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NonlinearMaximizer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RRMaximizer.cpp
    )

target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
