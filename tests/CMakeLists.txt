find_package( Catch2 3.8.0 CONFIG)
if( NOT Catch2_FOUND )
  message( STATUS "Could NOT Find Catch2 (Building v3.8.0 Locally)" )
  include( FetchContent )
  FetchContent_Declare(
    catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.8.0
  )
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(catch2)
else()
  message( STATUS "Found Catch2 VERSION ${Catch2_VERSION} DIR ${Catch2_DIR}" )
endif()

add_executable(mrchem-tests unit_tests.cpp)

add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)
add_subdirectory(solventeffect)

target_link_libraries(mrchem-tests
    PUBLIC
    mrchem
    Catch2::Catch2WithMain
  )

# Integration tests
add_subdirectory(h_el_field)
add_subdirectory(h2_scf_hf)
add_subdirectory(h2_pol_lda)
add_subdirectory(h2_mag_lda)
add_subdirectory(h2o_energy_blyp)
add_subdirectory(h2o_hirshfeld_lda)
add_subdirectory(h2_pol_cube)
add_subdirectory(li_scf_pbe0)
add_subdirectory(li_pol_lda)
add_subdirectory(hf_grad_lda)
add_subdirectory(hf_grad_blyp_surface_force)
add_subdirectory(cube_parser)
add_subdirectory(h2_scf_cube)
add_subdirectory(li_solv)
add_subdirectory(na_psp)
add_subdirectory(h_pb)
add_subdirectory(h_lpb)
add_subdirectory(he_zora_scf_lda)
add_subdirectory(he_azora_scf_lda)
add_subdirectory(cavity_input_parser)
add_subdirectory(h2_pol_solv)
