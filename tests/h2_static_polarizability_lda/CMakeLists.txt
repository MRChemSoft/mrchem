if(ENABLE_MPI)
    set(_h2_pol_lda_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_static_polarizability_LDA"
  LABELS "mrchem;H2;lda;scf;polarizability;short"
  COST 100
  LAUNCH_AGENT ${_h2_pol_lda_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
