if(ENABLE_MPI)
    set(_li_solv_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Li_solvent_effect"
  LABELS "mrchem;li_solv;solvent;scf;energy"
  COST 100
  LAUNCH_AGENT ${_li_solv_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
