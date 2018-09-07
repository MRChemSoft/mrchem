# Unit tests
if(ENABLE_UNIT_TESTS)
  # FIXME
  link_directories(${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

  # Prepare Catch "library"
  add_library(Catch INTERFACE)
  target_include_directories(Catch INTERFACE ${CMAKE_SOURCE_DIR}/external/catch)

  add_executable(mrchem-tests unit_tests.cpp)

  add_subdirectory(qmfunctions)
  add_subdirectory(qmoperators)

  # FIXME target mrchem must define its usage requirements
  target_include_directories(mrchem-tests
    PRIVATE
      ${CMAKE_BINARY_DIR}
      ${CMAKE_SOURCE_DIR}/src
      ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}
      ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/MRCPP
    )

  target_link_libraries(mrchem-tests
    PRIVATE
      mrchem
      Catch
    )
endif()

# Integration tests
if(ENABLE_TESTS)
  add_subdirectory(he)
  #add_subdirectory(li)
  #add_subdirectory(h2)
  #add_subdirectory(h2o)
  # This one is way too long...
  #add_subdirectory(c6h6)
endif()
