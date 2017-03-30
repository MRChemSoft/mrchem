include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwanalytic)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwcore)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwtrees)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwoperators)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwbuilders)

add_library(qmoperators STATIC 
    CoulombPotential.cpp
    DeltaOperator.cpp
    ExchangePotential.cpp
    FockOperator.cpp
    KineticOperator.cpp
    QMOperator.cpp
    QMPotential.cpp
    QMDerivative.cpp
    QMOperatorExp.cpp
    SpinOperator.cpp
    XCFunctional.cpp
    XCPotential.cpp
    XCOperator.cpp
)

target_link_libraries(qmoperators
    mrchem
    chemistry
    initial_guess
    qmfunctions
    mrcpp
    mwanalytic
    mwcore
    mwtrees
    mwoperators
    mwbuilders
    xcfun
)

install(TARGETS qmoperators
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

