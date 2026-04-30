# chemtensor doesn't support find_package yet
message(STATUS "Downloading and building chemtensor.")
include(FetchContent)
FetchContent_Declare(chemtensor
  QUIET
	GIT_REPOSITORY https://github.com/qc-tum/chemtensor.git
	GIT_TAG 18c269a5f1c4844c6d68188baa515aa4c9329611
)

# We save CMAKE_BUILD_TYPE, as we will set it to Release for externals
set(_build_type ${CMAKE_BUILD_TYPE})
set(_build_shared ${BUILD_SHARED_LIBS})
set(CMAKE_BUILD_TYPE "Debug")
set(BUILD_SHARED_LIBS OFF)

set(ENABLE_CHEMTENSOR_UNIT_TESTS OFF CACHE BOOL "" FORCE)
set(ENABLE_CHEMTENSOR_PYTHON_MODULE OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(chemtensor)

# Restore the settings
set(CMAKE_BUILD_TYPE ${_build_type})
set(BUILD_SHARED_LIBS ${_build_shared})