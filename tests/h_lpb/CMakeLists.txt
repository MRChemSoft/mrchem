if(ENABLE_MPI)
    set(_h_lpb_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H_linear_poisson_boltzmann"
  LABELS "mrchem;h_lpb;linear_poisson_boltzmann;scf;energy"
  COST 100
  LAUNCH_AGENT ${_h_lpb_launcher}
  )
