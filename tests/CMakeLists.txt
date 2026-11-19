option(ENABLE_TESTS "Enable test suite" ON)

if(ENABLE_TESTS)
option(CATCH2_FIND_BEHAVIOUR "Find Catch2 behavior" "default")
cpm_set_find_behaviour(${CATCH2_FIND_BEHAVIOUR})
CPMAddPackage(
  NAME Catch2
  VERSION 3
  GIT_TAG v3.16.0
  GITHUB_REPOSITORY catchorg/Catch2
  FIND_PACKAGE_ARGUMENTS "CONFIG"
  OPTIONS
  "BUILD_TESTING OFF"
  )

add_executable(mrchem-tests unit_tests.cpp)

# Suppress warnings from Eigen headers (-isystem instead of -I)
target_include_directories(mrchem-tests SYSTEM PRIVATE
  $<TARGET_PROPERTY:Eigen3::Eigen,INTERFACE_INCLUDE_DIRECTORIES>
  )

add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)
add_subdirectory(solventeffect)

target_link_libraries(mrchem-tests
    PUBLIC
    mrchem
    Catch2::Catch2WithMain
  )

# Integration tests
add_subdirectory(2h2o_far_energy_pbe)
# add_subdirectory(beh2_snqm_geomopt)
add_subdirectory(cavity_input_parser)
add_subdirectory(cube_parser)
add_subdirectory(h_linear_poisson_boltzmann)
add_subdirectory(h_poisson_boltzmann)
add_subdirectory(h_scf_finite_electric_field)
add_subdirectory(h2_magnetic_properties_lda)
add_subdirectory(h2_polarizability_cube)
add_subdirectory(h2_polarizability_solvent_effect)
add_subdirectory(h2_scf_cube)
add_subdirectory(h2_scf_hf)
add_subdirectory(h2_scf_libxc_b3lyp)
add_subdirectory(h2_scf_libxc_lda)
add_subdirectory(h2_scf_xcfun_b3lyp)
add_subdirectory(h2_scf_xcfun_lda)
add_subdirectory(h2_static_polarizability_lda)
# add_subdirectory(h2_znum)
add_subdirectory(h2o_energy_blyp)
add_subdirectory(h2o_hirshfeld_lda)
add_subdirectory(he_azora_scf_lda)
add_subdirectory(he_zora_scf_lda)
add_subdirectory(hf_grad_blyp_surface_force)
add_subdirectory(hf_molecular_gradient_lda)
add_subdirectory(li_scf_libxc_b3lyp)
add_subdirectory(li_scf_libxc_lda)
add_subdirectory(li_scf_xcfun_lda)
add_subdirectory(li_scf_xcfun_pbe0)
add_subdirectory(li_solvent_effect)
add_subdirectory(li_static_polarizability_lda)
add_subdirectory(li2h_imom_res)
add_subdirectory(li2h_imom_unres)
add_subdirectory(lih_core_hole)
add_subdirectory(na_pp_scf_hf)
add_subdirectory(ne_zora_scf_pbe)
add_subdirectory(o2_gto_qzvp_libxc_pbe)

endif()
