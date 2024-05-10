if(ENABLE_MPI)
    set(_h2o_core_hole_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2O_core_hole"
  LABELS "mrchem;h2o_core_hole;H2O_core_hole;dft;gga;energy"
  COST 100
  LAUNCH_AGENT ${_h2o_core_hole_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
