if(ENABLE_MPI)
    set(_h2_scf_svwn5_libxc "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_SCF_SVWN5_Libxc"
  LABELS "mrchem;h2_scf_svwn5_libxc;libxc;H2_SCF_SVWN5_Libxc;dft;scf;energy"
  COST 100
  LAUNCH_AGENT ${_h2_scf_svwn5_libxc_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
