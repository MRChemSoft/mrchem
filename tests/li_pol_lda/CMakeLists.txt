if(ENABLE_MPI)
    set(_li_pol_lda_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Li_static_polarizability_LDA"
  LABELS "li_pol_lda;mrchem;open_shell;polarizability;Li_static_polarizability_LDA"
  COST 100
  LAUNCH_AGENT ${_li_pol_lda_launcher}
  )
