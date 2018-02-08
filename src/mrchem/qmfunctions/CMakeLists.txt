include_directories (${CMAKE_SOURCE_DIR}/src/mrchem)

add_library(qmfunctions STATIC 
    QMFunction.cpp
    Density.cpp
    DensityProjector.cpp
    NuclearFunction.cpp
    Orbital.cpp
    OrbitalAdder.cpp
    OrbitalMultiplier.cpp
    OrbitalProjector.cpp
    OrbitalVector.cpp
)

target_link_libraries(qmfunctions
    mrchem
    chemistry
    initial_guess
    mrcpp
)

install(TARGETS qmfunctions
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

