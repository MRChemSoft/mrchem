if(ENABLE_MPI)
    set(_h2_tda_cube_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_tda_excited_states_CUBE"
  LABELS "H2_tda_excited_states_CUBE;excited_states;mrchem;h2_tda_cube"
  COST 100
  LAUNCH_AGENT ${_h2_tda_cube_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
