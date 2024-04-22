if(ENABLE_MPI)
    set(_h_pb_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H_poisson_boltzmann"
  LABELS "mrchem;h_pb;poisson_boltzmann;scf;energy"
  COST 100
  LAUNCH_AGENT ${_h_pb_launcher}
  )
