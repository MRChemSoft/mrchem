cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(MRChem_core
  VERSION 1.0.0
  LANGUAGES CXX
  )

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
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

add_subdirectory(analyticfunctions)
add_subdirectory(chemistry)
add_subdirectory(mrdft)
add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)
add_subdirectory(scf_solver)
add_subdirectory(utils)
add_subdirectory(initial_guess)

list(APPEND _public_headers
  SCFDriver.h
  mrchem.h
  mrenv.h
  parallel.h
  ${CONFIG_H_LOCATION}/config.h
  )
target_sources(mrchem
  PRIVATE
    SCFDriver.cpp
    parallel.cpp
    mrenv.cpp
  PUBLIC
    $<BUILD_INTERFACE:$<JOIN:${CMAKE_CURRENT_LIST_DIR}/,${_public_headers}>>
  )

target_link_libraries(mrchem
  PUBLIC
    Eigen3::Eigen
    MRCPP::mrcpp
    XCFun::xcfun
    getkw::getkw
  )

set_target_properties(mrchem
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    MACOSX_RPATH ON
    PUBLIC_HEADER "${_public_headers}"
  )

set(INSTALL_CMAKEDIR "share/cmake/MRChem")

install(
  TARGETS
    mrchem
  EXPORT
    MRChemTargets
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT lib
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT bin
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT lib
  PUBLIC_HEADER
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MRChem
    COMPONENT dev
  )

install(
  EXPORT
    MRChemTargets
  NAMESPACE
    "MRChem::"
  DESTINATION
    ${INSTALL_CMAKEDIR}
  COMPONENT
    dev
  )

add_executable(mrchem.x mrchem.cpp)

target_link_libraries(mrchem.x
  PRIVATE
    mrchem
  )

install(
  TARGETS
    mrchem.x
  EXPORT
    MRChemTargets-mrchemx
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT bin
  )

install(
  EXPORT
    MRChemTargets-mrchemx
  NAMESPACE
    "MRChem::"
  DESTINATION
    ${INSTALL_CMAKEDIR}
  COMPONENT
    bin
  )

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/MRChemConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
  )

configure_package_config_file(
  ${CMAKE_CURRENT_LIST_DIR}/MRChemConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/MRChemConfig.cmake
  INSTALL_DESTINATION ${INSTALL_CMAKEDIR}
  )

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/MRChemConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/MRChemConfigVersion.cmake
  DESTINATION
    ${INSTALL_CMAKEDIR}
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

macro(add_catch_test _name _labels)
  # _labels is not a list, it's a string... Transform it into a list
  set(labels)
  string(REPLACE ";" " " _labels "${_labels}")
  foreach(_label "${_labels}")
    list(APPEND labels ${_label})
  endforeach()
  unset(_labels)

  add_test(NAME ${_name}
           COMMAND ${PROJECT_BINARY_DIR}/tests/mrchem-tests [${_name}] --success --out ${PROJECT_BINARY_DIR}/tests/${_name}.log --durations yes
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  if(labels)
    set_tests_properties(${_name} PROPERTIES LABELS "${labels}")
  endif()
endmacro()


# Prepare Catch "library" for other executables
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/../external/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

enable_testing()
include(CTest)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../tests ${CMAKE_CURRENT_BINARY_DIR}/tests)
