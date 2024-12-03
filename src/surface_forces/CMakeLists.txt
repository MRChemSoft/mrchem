target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/SurfaceForce.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lebedev.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LebedevData.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/detail/lebedev_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/xcStress.cpp
    )
