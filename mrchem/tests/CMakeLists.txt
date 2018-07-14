# Prepare Catch "library"
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${MRChem_SOURCE_ROOT}/external/catch)

add_executable(mrchem-tests unit_tests.cpp)

add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)

target_link_libraries(mrchem-tests
  PRIVATE
    mrchem
    Catch
  )
