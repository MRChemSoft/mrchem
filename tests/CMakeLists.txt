# Unit tests
# FIXME
link_directories(${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

# Prepare Catch "library"
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CMAKE_SOURCE_DIR}/external/catch)

add_executable(mrchem-tests unit_tests.cpp)

add_subdirectory(qmfunctions)
#add_subdirectory(qmoperators)

target_link_libraries(mrchem-tests
  PRIVATE
    mrchem
    Catch
  )

# Integration tests
#add_subdirectory(h_el_field)
#add_subdirectory(h2_scf_hf)
#add_subdirectory(h2_pol_lda)
#add_subdirectory(h2o_energy_blyp)
#add_subdirectory(li_scf_pbe0)
#add_subdirectory(li_pol_lda)
