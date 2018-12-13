# Much disgusting
link_directories(${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

# <<< Build MRChem library >>>
add_library(mrchem
  parallel.cpp
  mrenv.cpp
  )

target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

# This command will become unnecessary, eventually
target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>
    $<BUILD_INTERFACE:${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP>
  )

list(APPEND _public_headers
  #SCFDriver.h
  mrchem.h
  mrenv.h
  parallel.h
  ${CMAKE_BINARY_DIR}/config.h
  )

add_dependencies(mrchem
  mrcpp_external
  xcfun_external
  )

add_subdirectory(chemistry)
add_subdirectory(analyticfunctions)
add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)
add_subdirectory(mrdft)
#add_subdirectory(scf_solver)
add_subdirectory(utils)

target_link_libraries(mrchem
  PUBLIC
    Eigen3::Eigen
    mrcpp
    xcfun
    getkw
  )

set_target_properties(mrchem
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    MACOSX_RPATH ON
    OUTPUT_NAME "mrchem"
    PUBLIC_HEADER "${_public_headers}"
  )

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

target_sources(mrchem.x
  PRIVATE
  #  SCFDriver.cpp
  PUBLIC
    $<BUILD_INTERFACE:$<JOIN:${CMAKE_CURRENT_LIST_DIR}/,${_public_headers}>>
  )

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

add_dependencies(mrchem.x
  mrcpp_external
  xcfun_external
  getkw_external
  )

target_link_libraries(mrchem.x
  PRIVATE
    mrchem
    mrcpp
    getkw
    xcfun
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
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT bin
  )

# <<< Build initial guess executables >>>
#add_subdirectory(initial_guess)

set(MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/mrchem.x)
# configure into CMakeFiles directory, because it must be in a different directory
# since file(COPY) does rename files
configure_file(mrchem.in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mrchem @ONLY)
# now copy the temporary into the final destination, setting the permissions
file(
  COPY ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mrchem
  DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
  FILE_PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
  )

install(
  PROGRAMS
    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/mrchem
  DESTINATION
    ${CMAKE_INSTALL_BINDIR}
  )
