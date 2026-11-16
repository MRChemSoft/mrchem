# HDF5 is optional: it is only needed to read GTO data written by DIRAC, Serenity or VeloxChem.
# With the default "none" it is neither searched for, fetched nor linked.
set(HDF5_FIND_BEHAVIOUR "none" CACHE STRING "Find HDF5 behavior (none, default, onlylocal, onlyfetch)")
set(ENABLE_HDF5 OFF)
set(MRCHEM_HDF5_LIBRARIES "")

if(HDF5_FIND_BEHAVIOUR STREQUAL "none")
  message(STATUS "HDF5 support disabled (enable with --hdf5=default)")
else()
  # HDF5 declares compatibility only within one minor version: an installed 2.1.x rejects
  # a request for 2, 2.0, 2.2 or a range like 2.0...<3. To accept any installed HDF5 2.x,
  # try each minor version in turn (a rejected request does not import any targets).
  if(HDF5_FIND_BEHAVIOUR STREQUAL "default" OR HDF5_FIND_BEHAVIOUR STREQUAL "onlylocal")
    foreach(_minor RANGE 9 0 -1)
      if(NOT HDF5_FOUND)
        find_package(HDF5 2.${_minor} CONFIG QUIET
          NO_CMAKE_PATH
          NO_CMAKE_PACKAGE_REGISTRY
          NO_CMAKE_SYSTEM_PACKAGE_REGISTRY
          )
      endif()
    endforeach()
  endif()

  if(HDF5_FOUND)
    message(STATUS "Found HDF5 ${HDF5_VERSION}: ${HDF5_DIR}")
  elseif(HDF5_FIND_BEHAVIOUR STREQUAL "onlylocal")
    message(FATAL_ERROR "--hdf5=onlylocal was given, but no installed HDF5 2.x was found")
  else()
    # HDF5 2.x itself requires CMake >= 3.26; MRChem only needs it when building HDF5 from source
    if(CMAKE_VERSION VERSION_LESS 3.26)
      message(FATAL_ERROR "Building HDF5 from source requires CMake >= 3.26 (found ${CMAKE_VERSION}). "
                          "Either install HDF5 2.x (e.g. 'conda install -c conda-forge hdf5'), "
                          "upgrade CMake, or configure without HDF5 (--hdf5=none).")
    endif()
    # Already searched above, so only fetch and build
    cpm_set_find_behaviour(onlyfetch)
    CPMAddPackage(
      NAME HDF5
      VERSION 2.2.0
      GITHUB_REPOSITORY HDFGroup/hdf5
      GIT_TAG 2.2.0
      OPTIONS
      "BUILD_TESTING OFF"
      "BUILD_SHARED_LIBS ON"
      "BUILD_STATIC_LIBS OFF"
      "HDF5_BUILD_TOOLS OFF"
      "HDF5_BUILD_EXAMPLES OFF"
      "HDF5_BUILD_HL_LIB OFF"
      "HDF5_BUILD_CPP_LIB ON"
      "HDF5_BUILD_FORTRAN OFF"
      "HDF5_BUILD_JAVA OFF"
      )
  endif()

  # Built from source -> plain target names; found installed -> hdf5:: namespace.
  # Alias so the namespaced names can always be used.
  if(TARGET hdf5-shared AND NOT TARGET hdf5::hdf5-shared)
    add_library(hdf5::hdf5-shared ALIAS hdf5-shared)
  endif()
  if(TARGET hdf5_cpp-shared AND NOT TARGET hdf5::hdf5_cpp-shared)
    add_library(hdf5::hdf5_cpp-shared ALIAS hdf5_cpp-shared)
  endif()

  if(NOT TARGET hdf5::hdf5-shared OR NOT TARGET hdf5::hdf5_cpp-shared)
    message(FATAL_ERROR "HDF5 was requested (--hdf5=${HDF5_FIND_BEHAVIOUR}) but the shared C and C++ libraries "
                        "(hdf5-shared, hdf5_cpp-shared) are not available")
  endif()

  set(ENABLE_HDF5 ON)
  set(MRCHEM_HDF5_LIBRARIES hdf5::hdf5-shared hdf5::hdf5_cpp-shared)
  message(STATUS "HDF5 support enabled")
endif()
