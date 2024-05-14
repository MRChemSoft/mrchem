if(ENABLE_MPI)
    set(_lih_core_hole_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "LiH_core_hole"
  LABELS "mrchem;lih_core_hole;LiH_core_hole;dft;gga;energy"
  COST 100
  LAUNCH_AGENT ${_lih_core_hole_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
