if(ENABLE_MPI)
    set(_li_scf_pbe0_libxc_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Li_SCF_PBE0_Libxc"
  LABELS "mrchem;li_scf_pbe0;Li_SCF_PBE0;dft;scf;energy;libxc"
  COST 100
  LAUNCH_AGENT ${_li_scf_pbe0_libxc_launcher}
  )
