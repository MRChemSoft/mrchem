if(ENABLE_MPI)
    set(_ne_scf_b3lyp_libxc_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Ne_SCF_B3LYP_Libxc"
  LABELS "mrchem;ne_scf_b3lyp_libxc;Ne_SCF_B3LYP_Libxc;energy;b3lyp;scf;ne;neon"
  COST 100
  LAUNCH_AGENT ${_ne_scf_b3lyp_libxc_launcher}
  )
