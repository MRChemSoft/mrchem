if(ENABLE_MPI)
    set(_h2_scf_libxc_b3lyp_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_SCF_Libxc_B3LYP"
  LABELS "mrchem;H2_SCF_Libxc_B3LYP;libxc;dft;scf;energy"
  COST 100
  LAUNCH_AGENT ${_h2_scf_libxc_b3lyp_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
