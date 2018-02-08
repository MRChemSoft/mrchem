include_directories (${CMAKE_SOURCE_DIR}/src/mrchem/chemistry)
include_directories (${CMAKE_SOURCE_DIR}/src/mrchem/properties)
include_directories (${CMAKE_SOURCE_DIR}/src/mrchem/qmfunctions)
#include_directories (${CMAKE_SOURCE_DIR}/src/mrchem/qmoperators)
#include_directories (${CMAKE_SOURCE_DIR}/src/mrchem/scf_solver)

add_subdirectory (initial_guess)
add_subdirectory (chemistry)
add_subdirectory (qmfunctions)
#add_subdirectory (qmoperators)
#add_subdirectory (src_solver)

#add_library(mrchem STATIC 
#  parallel.cpp
#)
#add_library(mrchem STATIC 
    #mrchem.in
    #mrchem.cpp
    #MREnv.cpp
    #SCFDriver.cpp
#)

#target_link_libraries(mrchem 
    #chemistry 
    #initial_guess
    #qmfunctions
    #qmoperators
    #scf_solver
    #mrcpp 
    #xcfun
    #getkw
#)

#add_dependencies(mrchem external-modules)

#add_executable(mrchem.x mrchem.cpp)
#add_dependencies(mrchem.x external-modules)

#target_link_libraries(mrchem.x 
    #mrchem 
    #chemistry
    #initial_guess
    #qmfunctions
    #qmoperators
    #scf_solver
    #mrcpp 
    #xcfun
    #getkw
    #${BLAS_LIBRARIES}
#)

#option (MRCHEM_DRIVER_DEBUG "Debug the driver (extmod)" OFF)

#set (MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/bin/mrchem.x)
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mrchem.in mrchem)

#install(TARGETS mrchem.x
	#RUNTIME DESTINATION bin
	#)

#install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/mrchem
	#DESTINATION bin
	#)
