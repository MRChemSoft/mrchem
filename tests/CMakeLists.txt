# Unit tests
# FIXME
link_directories(${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

# Prepare Catch "library"
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CMAKE_SOURCE_DIR}/external/catch)

add_executable(mrchem-tests unit_tests.cpp)

add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)

target_link_libraries(mrchem-tests
  PRIVATE
    mrchem
    Catch
  )

# Integration tests
add_subdirectory(li)
add_subdirectory(h2)
add_subdirectory(h2o)
