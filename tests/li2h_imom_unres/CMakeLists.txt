if(ENABLE_MPI)
    set(_li2h_imom_unres_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

if(ENABLE_XCFUN)
add_integration_test(
  NAME "Li2H_imom_unres"
  LABELS "mrchem;Li2H;scf;dft;gga;pbe;unrestricted;imom;energy;xcfun;medium"
  COST 100
  LAUNCH_AGENT ${_li2h_imom_unres_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
endif()
