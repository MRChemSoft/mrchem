if(ENABLE_MPI)
    set(_h2_scf_libxc_svwn5 "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_SCF_Libxc_SVWN5"
  LABELS "mrchem;H2;dft;scf;lda;svwn5;energy;short"
  COST 100
  LAUNCH_AGENT ${_h2_scf_libxc_svwn5_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
