if(ENABLE_MPI)
    set(_li_solv_pb_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Li_solvent_effect_poisson_boltzmann"
  LABELS "mrchem;li_solv_pb;solvent_pb;scf;energy"
  COST 100
  LAUNCH_AGENT ${_li_solv_pb_launcher}
  )
