include(DefaultPathsMRChem)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwanalytic)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwcore)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwtrees)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwoperators)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwbuilders)
include_directories (${CMAKE_SOURCE_DIR}/src/mrchem)
include_directories (${CMAKE_SOURCE_DIR}/src/mrchem/chemistry)
include_directories (${CMAKE_SOURCE_DIR}/src/mrchem/initial_guess)
include_directories (${CMAKE_SOURCE_DIR}/src/mrchem/qmfunctions)
include_directories (${CMAKE_SOURCE_DIR}/src/mrchem/qmoperators)
include_directories (${CMAKE_SOURCE_DIR}/src/mrchem/scf_solver)
include_directories (${CMAKE_SOURCE_DIR}/src/mrchem/properties)

link_directories (${CMAKE_INSTALL_PREFIX}/lib)

set (MRCHEM_EXECUTABLE ${CMAKE_BINARY_DIR}/pilot/mrchem-pilot.x)
set (MRCHEM_DRIVER_DEBUG ON)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mrchem.cpp)
    add_executable(mrchem-pilot.x mrchem.cpp)
    add_dependencies(mrchem-pilot.x external-modules)

    target_link_libraries(mrchem-pilot.x 
        mrchem 
        chemistry 
        initial_guess
        qmfunctions
        qmoperators
        scf_solver
        mrcpp 
        mwanalytic 
        mwcore
        mwtrees 
        mwoperators 
        mwbuilders
        getkw
        xcfun
        ${BLAS_LIBRARIES}
    )
    configure_file(${CMAKE_SOURCE_DIR}/src/mrchem/mrchem.in ${CMAKE_BINARY_DIR}/pilot/mrchem)
    execute_process(COMMAND chmod 755 ${CMAKE_BINARY_DIR}/pilot/mrchem OUTPUT_QUIET)
endif()

if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/mrchem.inp)
    file(COPY mrchem.inp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
