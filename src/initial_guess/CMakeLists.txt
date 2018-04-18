link_directories(${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

##################################
# Install Hydrogen initial guess #
##################################

# Install sources in lib
list(APPEND hydrogenguess_cpp
    hydrogen/hydrogen_guess.cpp
)
add_library(hydrogenguess STATIC ${hydrogenguess_cpp})
add_dependencies(hydrogenguess external-modules)
target_include_directories(hydrogenguess PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/chemistry
    ${CMAKE_SOURCE_DIR}/src/properties
    ${CMAKE_SOURCE_DIR}/src/analyticfunctions
    ${CMAKE_SOURCE_DIR}/src/qmfunctions
    ${CMAKE_SOURCE_DIR}/src/qmoperators
    ${CMAKE_SOURCE_DIR}/src/qmoperators/one_electron
    ${CMAKE_SOURCE_DIR}/src/qmoperators/two_electron
    ${CMAKE_SOURCE_DIR}/src/initial_guess/hydrogen
    )
target_link_libraries(hydrogenguess
    mrchem
    mrcpp
    )
install(TARGETS hydrogenguess
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/mrchem/initial_guess)

# Install headers in include
list(APPEND hydrogenguess_hpp
    hydrogen/hydrogen_guess.h
)
foreach(_header ${hydrogenguess_hpp})
   install(FILES ${_header}
           DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/mrchem/initial_guess/hydrogen)
endforeach()

# Install executable in bin
add_executable(hydrogen-guess.x hydrogen.cpp)
add_dependencies(hydrogen-guess.x external-modules)
target_include_directories(hydrogen-guess.x PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
    ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/chemistry
    ${CMAKE_SOURCE_DIR}/src/properties
    ${CMAKE_SOURCE_DIR}/src/analyticfunctions
    ${CMAKE_SOURCE_DIR}/src/qmfunctions
    ${CMAKE_SOURCE_DIR}/src/qmoperators
    ${CMAKE_SOURCE_DIR}/src/qmoperators/one_electron
    ${CMAKE_SOURCE_DIR}/src/qmoperators/two_electron
    ${CMAKE_SOURCE_DIR}/src/initial_guess/hydrogen
    )
target_link_libraries(hydrogen-guess.x
    hydrogenguess
    mrchem
    mrcpp
    getkw
    )
install(TARGETS hydrogen-guess.x
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
