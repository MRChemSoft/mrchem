if(ENABLE_MPI)
    set(_li2h_imom_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Li2H_imom"
  LABELS "mrchem;li2h_imom;Li2H_imom;dft;gga;energy"
  COST 100
  LAUNCH_AGENT ${_li2h_imom_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
