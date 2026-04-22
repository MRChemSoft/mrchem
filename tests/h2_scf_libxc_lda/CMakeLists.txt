if(ENABLE_MPI)
    set(_h2_scf_libxc_lda "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_SCF_Libxc_LDA"
  LABELS "mrchem;H2;libxc;dft;scf;lda;svwn5;energy;short"
  COST 100
  LAUNCH_AGENT ${_h2_scf_libxc_lda_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
