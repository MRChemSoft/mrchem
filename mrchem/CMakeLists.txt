cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(MRChem_core
  VERSION 1.0.0
  LANGUAGES CXX
  )

list(APPEND CMAKE_MODULE_PATH ${MRChem_SOURCE_ROOT}/cmake/custom)

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

add_subdirectory(src)

option(MRCHEM_DRIVER_DEBUG "Debug the driver (extmod)" OFF)
set(MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/mrchem.x)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mrchem.in ${CMAKE_CURRENT_BINARY_DIR}/mrchem @ONLY)
install(
  PROGRAMS
    ${CMAKE_CURRENT_BINARY_DIR}/mrchem
  DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}
  )

if(ENABLE_UNIT_TESTS)
  include(tests)
  enable_testing()
  include(CTest)
  add_subdirectory(tests) # This must come last!!
endif()
