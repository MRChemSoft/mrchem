if(ENABLE_MPI)
    set(_he_zora_scf_lda_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "He_ZORA_SCF_LDA"
  LABELS "mrchem;He;dft;lda;scf;energy;zora;short"
  COST 100
  LAUNCH_AGENT ${_he_zora_scf_lda_launcher}
  )
