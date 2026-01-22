if(ENABLE_MPI)
    set(_he_scf_lda_libxc_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "He_SCF_LDA_LIBXC"
  LABELS "mrchem;he_scf_lda_libxc;He_SCF_LDA_LIBXC;energy;lda;scf;libxc;he;helium"
  COST 100
  LAUNCH_AGENT ${_he_scf_lda_libxc_launcher}
  )
