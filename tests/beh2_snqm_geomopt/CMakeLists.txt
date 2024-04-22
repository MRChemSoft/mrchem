if(ENABLE_MPI)
    set(_beh2_snqm_geomopt_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "BeH2_SNQM_geometry_optimization"
  LABELS "mrchem;beh2;snqm;geometry_optimization;scf;lda"
  COST 100
  LAUNCH_AGENT ${_beh2_snqm_geomopt_launcher}
  )
