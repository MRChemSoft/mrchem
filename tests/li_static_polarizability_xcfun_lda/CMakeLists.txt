if(ENABLE_MPI)
    set(_li_pol_lda_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Li_static_polarizability_XCFun_LDA"
  LABELS "mrchem;Li;scf;lda;dft;unrestricted;polarizability;xcfun;medium"
  COST 100
  LAUNCH_AGENT ${_li_pol_lda_launcher}
  )
