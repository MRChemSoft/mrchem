if(ENABLE_MPI)
  set(_n_scf_pbe0_libxc_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "N_SCF_PBE0_Libxc"
  LABELS "mrchem;n;scf;dft;pbe0;energy;libxc"
  COST 150
  LAUNCH_AGENT ${_n_scf_pbe0_libxc_launcher}
)
