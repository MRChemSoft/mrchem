# Unit tests
# Prepare "Catch" library for other executables, if not already available
if(NOT TARGET Catch)
  add_library(Catch INTERFACE)
  target_include_directories(Catch INTERFACE ${PROJECT_SOURCE_DIR}/external/catch)
endif()

add_executable(mrchem-tests unit_tests.cpp)

add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)
add_subdirectory(solventeffect)

target_link_libraries(mrchem-tests
  PRIVATE
    mrchem
    Catch
  )

# Integration tests
add_subdirectory(h_el_field)
add_subdirectory(h2_scf_hf)
add_subdirectory(h2_pol_lda)
add_subdirectory(h2o_energy_blyp)
add_subdirectory(li_scf_pbe0)
add_subdirectory(li_pol_lda)
add_subdirectory(hf_grad_lda)
add_subdirectory(cube_parser)
add_subdirectory(h2_scf_cube)
add_subdirectory(li_solv)
add_subdirectory(he_zora_scf_lda)
