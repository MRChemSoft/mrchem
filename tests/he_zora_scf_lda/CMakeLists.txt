if(ENABLE_MPI)
    set(_he_zora_scf_lda_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "He_ZORA_SCF_LDA"
  LABELS "mrchem;he_zora_scf_lda;He_ZORA_SCF_LDA;energy;lda;scf;zora;he;helium;short"
  COST 100
  LAUNCH_AGENT ${_he_zora_scf_lda_launcher}
  )
