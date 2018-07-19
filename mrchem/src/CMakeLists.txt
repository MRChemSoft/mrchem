# <<< Build MRChem library >>>
add_library(mrchem
  parallel.cpp
  mrenv.cpp
  )

target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${MRChem_BINARY_ROOT}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/mrchem>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

list(APPEND _public_headers
  SCFDriver.h
  mrchem.h
  mrenv.h
  parallel.h
  ${MRChem_BINARY_ROOT}/config.h
  )

add_subdirectory(mrchem)

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
    OUTPUT_NAME "mrchem"
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

# <<< Build MRChem executables >>>
add_executable(mrchem.x mrchem.cpp)

target_sources(mrchem.x
  PRIVATE
    SCFDriver.cpp
  PUBLIC
    $<BUILD_INTERFACE:$<JOIN:${CMAKE_CURRENT_LIST_DIR}/,${_public_headers}>>
  )

target_link_libraries(mrchem.x
  PRIVATE
    mrchem
  )

# RPATH fixing
file(RELATIVE_PATH _rel ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR} ${CMAKE_INSTALL_PREFIX})
if(APPLE)
  set(_rpath "@loader_path/${_rel}")
else()
  set(_rpath "\$ORIGIN/${_rel}")
endif()
file(TO_NATIVE_PATH "${_rpath}/${CMAKE_INSTALL_LIBDIR}" MRChem_RPATH)
set_target_properties(mrchem.x
  PROPERTIES
    MACOSX_RPATH ON
    SKIP_BUILD_RPATH OFF
    BUILD_WITH_INSTALL_RPATH OFF
    INSTALL_RPATH "${MRChem_RPATH}"
    INSTALL_RPATH_USE_LINK_PATH ON
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

# <<< Build initial guess executables >>>
add_subdirectory(initial_guess)

# <<< Generate CMake configuration module >>>
#include(CMakePackageConfigHelpers)
#write_basic_package_version_file(
#  ${CMAKE_CURRENT_BINARY_DIR}/MRChemConfigVersion.cmake
#  VERSION ${PROJECT_VERSION}
#  COMPATIBILITY SameMajorVersion
#  )
#
#configure_package_config_file(
#  ${CMAKE_CURRENT_LIST_DIR}/MRChemConfig.cmake.in
#  ${CMAKE_CURRENT_BINARY_DIR}/MRChemConfig.cmake
#  INSTALL_DESTINATION ${INSTALL_CMAKEDIR}
#  )
#
#install(
#  FILES
#    ${CMAKE_CURRENT_BINARY_DIR}/MRChemConfig.cmake
#    ${CMAKE_CURRENT_BINARY_DIR}/MRChemConfigVersion.cmake
#  DESTINATION
#    ${INSTALL_CMAKEDIR}
#  )
