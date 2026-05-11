if(ENABLE_MPI)
    set(_li_scf_lda_xcfun_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Li_SCF_XCFun_LDA"
  LABELS "mrchem;Li;dft;lda;scf;unrestricted;xcfun;energy;short"
  COST 100
  LAUNCH_AGENT ${_li_scf_xcfun_svwn5_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
