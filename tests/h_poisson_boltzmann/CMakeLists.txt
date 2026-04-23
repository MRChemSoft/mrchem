if(ENABLE_MPI)
    set(_h_pb_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H_poisson_boltzmann"
  LABELS "mrchem;H;scf;pbe0;gga;pcm;poisson_boltzmann;energy;xcfun;short"
  COST 100
  LAUNCH_AGENT ${_h_pb_launcher}
  )
