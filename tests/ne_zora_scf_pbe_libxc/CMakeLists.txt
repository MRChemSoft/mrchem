if(ENABLE_MPI)
    set(_ne_zora_scf_pbe_libxc_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Ne_ZORA_SCF_PBE_Libxc"
  LABELS "mrchem;ne_zora_scf_pbe_libxc;Ne_ZORA_SCF_PBE_Libxc;energy;pbe;scf;zora;ne;neon"
  COST 100
  LAUNCH_AGENT ${_ne_zora_scf_pbe_libxc_launcher}
  )
