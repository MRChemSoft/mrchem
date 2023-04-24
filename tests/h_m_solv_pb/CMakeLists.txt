if(ENABLE_MPI)
    set(_h_m_solv_pb_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H_m_solvent_effect_poisson_boltzmann"
  LABELS "mrchem;h_m_solv_pb;solvent_pb;scf;energy"
  COST 100
  LAUNCH_AGENT ${_h_m_solv_pb_launcher}
  )
