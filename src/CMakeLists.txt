cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(MRChem_core LANGUAGES CXX)

include(GNUInstallDirs)

find_package(MRCPP CONFIG REQUIRED)
get_property(_loc TARGET MRCPP::mrcpp PROPERTY LOCATION)
message(STATUS "Using MRCPP: ${_loc} (version ${MRCPP_VERSION})")
# FIXME Uncomment these
#find_package(Getkw CONFIG REQUIRED)
#get_property(_loc TARGET Getkw::getkw PROPERTY LOCATION)
#message(STATUS "Using Getkw: ${_loc} (version ${Getkw_VERSION})")
#find_package(XCFun CONFIG REQUIRED)
#get_property(_loc TARGET XCFun::xcfun PROPERTY LOCATION)
#message(STATUS "Using XCFun: ${_loc} (version ${XCFun_VERSION})")

add_library(mrchem 
    SCFDriver.cpp
    parallel.cpp
    mrenv.cpp
    )
add_subdirectory(chemistry)
add_subdirectory(analyticfunctions)
add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)
add_subdirectory(mrdft)
add_subdirectory(scf_solver)
add_subdirectory(initial_guess)
add_subdirectory(utils)

target_include_directories(mrchem
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/chemistry
    ${CMAKE_CURRENT_SOURCE_DIR}/properties
    ${CMAKE_CURRENT_SOURCE_DIR}/analyticfunctions
    ${CMAKE_CURRENT_SOURCE_DIR}/qmfunctions
    ${CMAKE_CURRENT_SOURCE_DIR}/qmoperators
    ${CMAKE_CURRENT_SOURCE_DIR}/qmoperators/one_electron
    ${CMAKE_CURRENT_SOURCE_DIR}/qmoperators/two_electron
    ${CMAKE_CURRENT_SOURCE_DIR}/scf_solver
    ${CMAKE_CURRENT_SOURCE_DIR}/initial_guess
    )
target_include_directories(mrchem
  SYSTEM
  PUBLIC
    # FIXME Remove as soon as moving to Eigen 3.3
    ${EIGEN3_INCLUDE_DIR}
    # FIXME Remove as soon as updating to the latest and greatest XCFun
    ${XCFun_INCLUDEDIR}
    # FIXME Remove as soon as updating to the latest and greatest Getkw
    ${Getkw_INCLUDEDIR}
    )
target_link_libraries(mrchem
  PUBLIC
    MRCPP::mrcpp
    )
install(TARGETS mrchem ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/mrchem)

add_executable(mrchem.x
    mrchem.cpp
    )
target_link_libraries(mrchem.x
    mrchem
    hydrogenguess
    gtoguess
    # FIXME
    #XCFun::xcfun
    ${XCFun_LIBRARIES}
    # FIXME
    #Getkw::getkw
    ${Getkw_LIBRARIES}
    )
install(TARGETS mrchem.x RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

option(MRCHEM_DRIVER_DEBUG "Debug the driver (extmod)" OFF)
set(MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/mrchem.x)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mrchem.in ${CMAKE_CURRENT_BINARY_DIR}/mrchem)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/mrchem DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
