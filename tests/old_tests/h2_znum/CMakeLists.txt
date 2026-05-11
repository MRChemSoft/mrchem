if(ENABLE_MPI)
    set(_h2_znum_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_Z_number"
  LABELS "mrchem;h2_znum;H2_Z_number;energy;atomic_number;"
  COST 100
  LAUNCH_AGENT ${_h2_znum_launcher}
  )
