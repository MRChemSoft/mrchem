add_library(initial_guess STATIC 
            AOBasis.cpp
            AOContraction.cpp
            Intgrl.cpp
            OrbitalExp.cpp
            )

target_link_libraries(initial_guess mrcpp)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lsdalton.py.in lsdalton.py)

install(FILES 
        ${CMAKE_CURRENT_BINARY_DIR}/lsdalton.py
        periodictable.py
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}
        )
