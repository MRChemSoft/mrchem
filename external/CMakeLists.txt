include(ExternalProject)
add_custom_target(external-modules)

add_custom_target(git-submodule-init
    COMMAND git submodule init
    COMMAND git submodule update
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

ExternalProject_Add(submodule-mrcpp
    DOWNLOAD_COMMAND ""
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/mrcpp
    CMAKE_ARGS  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                -DCMAKE_BUILD_TYPE=Release
                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -DENABLE_OPENMP=${ENABLE_OPENMP}
                -DENABLE_MPI=${ENABLE_MPI}
                -DENABLE_TESTS=FALSE
                -DENABLE_EXAMPLES=FALSE
    CMAKE_CACHE_ARGS
                -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    INSTALL_COMMAND
                DESTDIR=${CMAKE_BINARY_DIR}/stage ${CMAKE_MAKE_PROGRAM} install
    )
add_dependencies(submodule-mrcpp git-submodule-init)
add_dependencies(external-modules submodule-mrcpp)

enable_language(C)
ExternalProject_Add(submodule-libgetkw
    DOWNLOAD_COMMAND ""
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libgetkw
    CMAKE_ARGS  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                -DCMAKE_BUILD_TYPE=Release
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    CMAKE_CACHE_ARGS
                -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    INSTALL_COMMAND
                DESTDIR=${CMAKE_BINARY_DIR}/stage ${CMAKE_MAKE_PROGRAM} install
    )
add_dependencies(submodule-libgetkw git-submodule-init)
add_dependencies(external-modules submodule-libgetkw)

option(ENABLE_XCFUN "Use XCFun library" ON)
if (ENABLE_XCFUN)
    ExternalProject_Add(submodule-xcfun
        DOWNLOAD_COMMAND ""
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/xcfun
        CMAKE_ARGS  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                    -DCMAKE_BUILD_TYPE=Release
                    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        CMAKE_CACHE_ARGS
                    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
        INSTALL_COMMAND
                    DESTDIR=${CMAKE_BINARY_DIR}/stage ${CMAKE_MAKE_PROGRAM} install
        )
    add_dependencies(submodule-xcfun git-submodule-init)
    add_dependencies(external-modules submodule-xcfun)
endif()
