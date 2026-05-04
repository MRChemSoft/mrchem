if(ENABLE_MPI)
    set(_ne_zora_scf_pbe_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Ne_ZORA_SCF_PBE"
  LABELS "mrchem;Ne;dft;pbe;gga;zora;energy;xcfun;medium"
  COST 100
  LAUNCH_AGENT ${_ne_zora_scf_pbe_launcher}
  )
