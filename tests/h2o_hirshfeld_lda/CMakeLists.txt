if(ENABLE_MPI)
    set(_h2o_hirshfeld_lda_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2O_Hirshfeld_charges_LDA"
  LABELS "mrchem;H2O;dft;lda;hirshfeld;charges;xcfun;short"
  COST 100
  LAUNCH_AGENT ${_h2o_hirshfeld_lda_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
