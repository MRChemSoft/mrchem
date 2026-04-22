if(ENABLE_MPI)
    set(_h2_scf_xcfun_lda "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_SCF_XCFun_LDA"
  LABELS "mrchem;H2;xcfun;dft;scf;lda;svwn5;energy;short"
  COST 100
  LAUNCH_AGENT ${_h2_scf_xcfun_xcfun_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
