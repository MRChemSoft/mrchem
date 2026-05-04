if(ENABLE_MPI)
    set(_2h2o_far_energy_pbe_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "2H2O_far_energy_pbe"
  LABELS "mrchem;H2O;dft;gga;pbe;xcfun;energy;medium"
  COST 100
  LAUNCH_AGENT ${_2h2o_far_energy_pbe_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
