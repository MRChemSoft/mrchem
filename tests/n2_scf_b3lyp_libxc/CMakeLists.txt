if(ENABLE_MPI)
  set(_n2_scf_b3lyp_libxc_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "N2_SCF_B3LYP_Libxc"
  LABELS "mrchem;n2;scf;dft;b3lyp;energy;libxc"
  COST 150
  LAUNCH_AGENT ${_n2_scf_b3lyp_libxc_launcher}
)
