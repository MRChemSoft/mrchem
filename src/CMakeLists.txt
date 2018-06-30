cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(MRChem_core LANGUAGES CXX)

include(GNUInstallDirs)

find_package(Eigen3 3.3 CONFIG REQUIRED)
message(STATUS "Using Eigen3: ${EIGEN3_ROOT_DIR} (version ${Eigen3_VERSION})")

find_package(MRCPP CONFIG REQUIRED)
get_property(_loc TARGET MRCPP::mrcpp PROPERTY LOCATION)
message(STATUS "Using MRCPP: ${_loc} (version ${MRCPP_VERSION})")

find_package(getkw CONFIG REQUIRED)
get_property(_loc TARGET getkw::getkw PROPERTY LOCATION)
message(STATUS "Using getkw: ${_loc} (version ${getkw_VERSION})")

find_package(XCFun CONFIG REQUIRED)
get_property(_loc TARGET XCFun::xcfun PROPERTY LOCATION)
message(STATUS "Using XCFun: ${_loc} (version ${XCFun_VERSION})")

add_library(mrchem "")

target_sources(mrchem
  PRIVATE
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
    ${CONFIG_H_LOCATION}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/chemistry
    ${CMAKE_CURRENT_SOURCE_DIR}/properties
    ${CMAKE_CURRENT_SOURCE_DIR}/analyticfunctions
    ${CMAKE_CURRENT_SOURCE_DIR}/qmfunctions
    ${CMAKE_CURRENT_SOURCE_DIR}/qmoperators
    ${CMAKE_CURRENT_SOURCE_DIR}/qmoperators/one_electron
    ${CMAKE_CURRENT_SOURCE_DIR}/qmoperators/two_electron
    ${CMAKE_CURRENT_SOURCE_DIR}/scf_solver
    ${CMAKE_CURRENT_SOURCE_DIR}/utils
  )
target_link_libraries(mrchem
  PUBLIC
    MRCPP::mrcpp
    Eigen3::Eigen
    getkw::getkw
    XCFun::xcfun
  )
install(
  TARGETS
    mrchem
  ARCHIVE DESTINATION
    ${CMAKE_INSTALL_LIBDIR}
  )

add_executable(mrchem.x
    mrchem.cpp
  )
target_link_libraries(mrchem.x
  PRIVATE
    mrchem
  )
install(
  TARGETS
    mrchem.x
  RUNTIME DESTINATION
    ${CMAKE_INSTALL_BINDIR}
  )

option(MRCHEM_DRIVER_DEBUG "Debug the driver (extmod)" OFF)
set(MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/mrchem.x)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mrchem.in ${CMAKE_CURRENT_BINARY_DIR}/mrchem @ONLY)
install(
  PROGRAMS
    ${CMAKE_CURRENT_BINARY_DIR}/mrchem
  DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}
  )
