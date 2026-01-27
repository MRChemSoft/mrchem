if(ENABLE_MPI)
    set(_he_scf_lda_libxc_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "He_SCF_LDA_LIBXC"
  LABELS "mrchem;energy;dft;scf;libxc"
  COST 100
  LAUNCH_AGENT ${_he_scf_lda_libxc_launcher}
  )
