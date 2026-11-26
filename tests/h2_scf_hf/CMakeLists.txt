if(ENABLE_MPI)
    set(_h2_scf_hf_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

if(ENABLE_XCFUN)
add_integration_test(
  NAME "H2_SCF_HF"
  LABELS "mrchem;H2;scf;HF;energy;xcfun;short"
  COST 100
  LAUNCH_AGENT ${_h2_scf_hf_launcher}
  )
endif()
