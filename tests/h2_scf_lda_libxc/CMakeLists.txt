if(ENABLE_MPI)
    set(_h2_scf_lda_libxc "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_SCF_LDA_Libxc"
  LABELS "mrchem;h2_scf_lda_libxc;libxc;H2_SCF_LDA_Libxc;dft;scf;energy"
  COST 100
  LAUNCH_AGENT ${_h2_scf_lda_libxc_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
