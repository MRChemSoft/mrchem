if(ENABLE_MPI)
    set(_h2_scf_b3lyp_libxc_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_SCF_B3LYP_Libxc"
  LABELS "mrchem;energy;dft;scf;libxc"
  COST 100
  LAUNCH_AGENT ${_h2_scf_b3lyp_libxc_launcher}
  )
