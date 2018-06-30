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

target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${CONFIG_H_LOCATION}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/properties>
  )

add_subdirectory(analyticfunctions)
add_subdirectory(chemistry)
add_subdirectory(mrdft)
add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)
add_subdirectory(scf_solver)
add_subdirectory(utils)
add_subdirectory(initial_guess)

target_sources(mrchem
  PRIVATE
    SCFDriver.cpp
    parallel.cpp
    mrenv.cpp
  )

target_link_libraries(mrchem
  PUBLIC
    Eigen3::Eigen
    MRCPP::mrcpp
    XCFun::xcfun
    getkw::getkw
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
