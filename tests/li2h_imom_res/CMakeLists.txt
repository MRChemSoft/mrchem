if(ENABLE_MPI)
    set(_li2h_imom_res_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Li2H_imom_res"
  LABELS "mrchem;li2h_imom_res;Li2H_imom_res;dft;gga;energy"
  COST 100
  LAUNCH_AGENT ${_li2h_imom_res_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
