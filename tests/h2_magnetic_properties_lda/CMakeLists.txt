if(ENABLE_MPI)
    set(_h2_mag_lda_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_magnetic_properties_LDA"
  LABELS "mrchem;H2;magnetizability;scf;nmr;lda;xcfun;short"
  COST 100
  LAUNCH_AGENT ${_h2_mag_lda_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
