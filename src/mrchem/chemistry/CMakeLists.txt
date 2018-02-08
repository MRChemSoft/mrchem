add_library(chemistry STATIC 
    Molecule.cpp
    PeriodicTable.cpp
)

target_link_libraries(chemistry
    #mrchem
    #qmfunctions
    #scf_solver
    mrcpp
)

install(TARGETS chemistry
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

