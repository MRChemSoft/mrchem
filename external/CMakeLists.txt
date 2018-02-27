include(ExternalProject)
add_custom_target(external-modules)

add_custom_target(git-submodule-init
    COMMAND git submodule init
    COMMAND git submodule update
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

enable_language(C)
ExternalProject_Add(submodule-libgetkw
    DOWNLOAD_COMMAND git submodule update
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}
    DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libgetkw
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libgetkw
    INSTALL_DIR ${PROJECT_BINARY_DIR}/external
    CMAKE_ARGS  -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
                -DCMAKE_BUILD_TYPE=Release
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    )
add_dependencies(submodule-libgetkw git-submodule-init)
add_dependencies(external-modules submodule-libgetkw)

option(ENABLE_XCFUN "Use XCFun library" ON)
if (ENABLE_XCFUN)
    ExternalProject_Add(submodule-xcfun
        DOWNLOAD_COMMAND git submodule update
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}
        DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/xcfun
        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/xcfun
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        CMAKE_ARGS  -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
                    -DCMAKE_BUILD_TYPE=Release
                    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        )
    add_dependencies(submodule-xcfun git-submodule-init)
    add_dependencies(external-modules submodule-xcfun)
endif()

include_directories(${PROJECT_SOURCE_DIR}/external/catch)
