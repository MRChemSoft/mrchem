if(ENABLE_MPI)
    set(_he_azora_scf_lda_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "He_AZORA_SCF_LDA"
  LABELS "mrchem;He;lda;dft;scf;zora;azora;energy;xcfun;short"
  COST 100
  LAUNCH_AGENT ${_he_azora_scf_lda_launcher}
  )
