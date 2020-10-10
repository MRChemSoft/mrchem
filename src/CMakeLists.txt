# <<< Build MRChem library >>>
set(BUILD_SHARED_LIBS OFF)
add_library(mrchem
  driver.cpp
  parallel.cpp
  mrenv.cpp
  )

add_dependencies(mrchem mrchem-info)

target_compile_definitions(mrchem
  PUBLIC
    $<$<AND:$<TARGET_EXISTS:OpenMP::OpenMP_CXX>,$<BOOL:${ENABLE_OPENMP}>>:MRCHEM_HAS_OMP>
    $<$<AND:$<TARGET_EXISTS:MPI::MPI_CXX>,$<BOOL:${ENABLE_MPI}>>:MRCHEM_HAS_MPI>
  )

target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

list(APPEND _public_headers
  ${CMAKE_CURRENT_SOURCE_DIR}/driver.h
  ${CMAKE_CURRENT_SOURCE_DIR}/mrchem.h
  ${CMAKE_CURRENT_SOURCE_DIR}/mrenv.h
  ${CMAKE_CURRENT_SOURCE_DIR}/parallel.h
  ${PROJECT_BINARY_DIR}/config.h
  ${PROJECT_BINARY_DIR}/version.h
  )

# RPATH fixing
file(RELATIVE_PATH _rel ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR} ${CMAKE_INSTALL_PREFIX})
if(APPLE)
  set(_rpath "@loader_path/${_rel}")
else()
  set(_rpath "\$ORIGIN/${_rel}")
endif()
file(TO_NATIVE_PATH "${_rpath}/${CMAKE_INSTALL_LIBDIR}" MRChem_RPATH)

add_subdirectory(chemistry)
add_subdirectory(analyticfunctions)
add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)
add_subdirectory(mrdft)
add_subdirectory(scf_solver)
add_subdirectory(utils)

target_link_libraries(mrchem
  PRIVATE
    Eigen3::Eigen
  PUBLIC
    XCFun::xcfun
    MRCPP::mrcpp
    nlohmann_json::nlohmann_json
    $<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>
    $<$<TARGET_EXISTS:MPI::MPI_CXX>:MPI::MPI_CXX>
    # NOTE the above should become this for consistency, but this requires #155 in MRCPP to be merged first
    #$<$<AND:$<TARGET_EXISTS:OpenMP::OpenMP_CXX>,$<BOOL:${ENABLE_OPENMP}>>:OpenMP::OpenMP_CXX>
    #$<$<AND:$<TARGET_EXISTS:MPI::MPI_CXX>,$<BOOL:${ENABLE_MPI}>>:MPI::MPI_CXX>
  )

set_target_properties(mrchem
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    MACOSX_RPATH ON
    OUTPUT_NAME "mrchem"
    PUBLIC_HEADER "${_public_headers}"
  )

if(APPLE)
  set_target_properties(mrchem
    PROPERTIES
      LINK_FLAGS "-undefined dynamic_lookup"
    )
endif()

install(
  TARGETS
    mrchem
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

# <<< Build MRChem executables >>>
add_executable(mrchem.x mrchem.cpp)

target_link_libraries(mrchem.x
  PRIVATE
    mrchem
  )

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
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT bin
  )

add_subdirectory(initial_guess)
