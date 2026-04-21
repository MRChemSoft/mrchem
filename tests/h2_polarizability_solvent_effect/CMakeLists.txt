if(ENABLE_MPI)
    set(_h2_pol_solv_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_polarizability_solvent_effect"
  LABELS "mrchem;H2;solvent;scf;HF;pcm;energy;polarizability;xcfun;medium"
  COST 150
  LAUNCH_AGENT ${_h2_pol_solv_launcher}
  )
