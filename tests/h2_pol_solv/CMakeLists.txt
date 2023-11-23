if(ENABLE_MPI)
    set(_h2_pol_solv_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_polarizability_solvent_effect"
  LABELS "mrchem;h2_pol_solv;solvent;scf;energy;polarizability"
  COST 150
  LAUNCH_AGENT ${_h2_pol_solv_launcher}
  )
