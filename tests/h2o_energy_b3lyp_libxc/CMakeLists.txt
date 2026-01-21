if(ENABLE_MPI)
    set(_h2o_energy_b3lyp_libxc_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2O_energy_B3LYP_libxc"
  LABELS "mrchem;h2o_energy_b3lyp_libxc;H2O_energy_B3LYP_libxc;dft;gga;energy"
  COST 100
  LAUNCH_AGENT ${_h2o_energy_b3lyp_libxc_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
