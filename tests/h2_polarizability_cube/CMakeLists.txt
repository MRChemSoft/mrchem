if(ENABLE_MPI)
    set(_h2_pol_cube_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_polarizability_CUBE"
  LABELS "mrchem;H2;lda;polarizability;cube;xcfun;medium"
  COST 100
  LAUNCH_AGENT ${_h2_pol_cube_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
