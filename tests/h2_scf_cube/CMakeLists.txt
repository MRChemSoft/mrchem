if(ENABLE_MPI)
    set(_h2_scf_cube_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_SCF_CUBE"
  LABELS "mrchem;h2_scf_cube;cube_guess;energy;mrchem;scf"
  COST 100
  LAUNCH_AGENT ${_h2_scf_cube_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
