include_directories (${CMAKE_INSTALL_PREFIX}/include) 

include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwanalytic)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwcore)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwtrees)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwoperators)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwbuilders)
include_directories (chemistry)
include_directories (initial_guess)
include_directories (qmfunctions)
include_directories (qmoperators)
include_directories (scf_solver)
include_directories (properties)

add_subdirectory (chemistry)
add_subdirectory (initial_guess)
add_subdirectory (qmfunctions)
add_subdirectory (qmoperators)
add_subdirectory (scf_solver)

link_directories (${CMAKE_INSTALL_PREFIX}/lib)

add_library(mrchem STATIC 
    mrchem.in
    mrchem.cpp
    MREnv.cpp
    SCFDriver.cpp
)

target_link_libraries(mrchem 
    chemistry 
    initial_guess
    qmfunctions
    qmoperators
    scf_solver
    mrcpp 
    mwanalytic
    mwtrees
    mwoperators
    mwbuilders
    getkw
    xcfun
)

add_dependencies(mrchem external-modules)

add_executable(mrchem.x mrchem.cpp)
add_dependencies(mrchem.x external-modules)

target_link_libraries(mrchem.x 
    mrchem 
    chemistry
    initial_guess
    qmfunctions
    qmoperators
    scf_solver
    mrcpp 
    mwanalytic
    mwtrees
    mwoperators
    mwbuilders
    getkw
    xcfun
    ${BLAS_LIBRARIES}
)

option (MRCHEM_DRIVER_DEBUG "Debug the driver (extmod)" OFF)

set (MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/bin/mrchem.x)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mrchem.in mrchem)

install(TARGETS mrchem.x
	RUNTIME DESTINATION bin
	)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/mrchem
	DESTINATION bin
	)
