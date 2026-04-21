if(ENABLE_MPI)
    set(_li_b3lyp_scf_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Li_SCF_Libxc_B3LYP"
  LABELS "mrchem;libxc;Li_SCF_Libxc_B3LYP;dft;scf;energy"
  COST 100
  LAUNCH_AGENT ${_li_scf_libxc_b3lyp_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
