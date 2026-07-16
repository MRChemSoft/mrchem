target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/XCLib.cpp
    )

if(ENABLE_LIBXC)
target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Libxc_func_alias.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Libxc.cpp
    )
endif()

if(ENABLE_XCFUN)
target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/XCFun.cpp
    )
endif()
