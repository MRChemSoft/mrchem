if(ENABLE_MPI)
    set(_h2_scf_b3lyp_libxc_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_SCF_B3LYP_Libxc"
  LABELS "mrchem;h2_scf_b3lyp_libxc;H2_SCF_B3LYP_Libxc;energy;b3lyp;scf;h2;hydrogen;libxc"
  COST 100
  LAUNCH_AGENT ${_h2_scf_b3lyp_libxc_launcher}
  )
