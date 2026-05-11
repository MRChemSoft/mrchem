if(ENABLE_MPI)
    set(_o2_gto_QZVP_libxc_pbe_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "O2_gto_QZVP_libxc_pbe"
  LABELS "mrchem;O2;dft;gga;pbe;gto;libxc;short"
  COST 100
  LAUNCH_AGENT ${_o2_gto_QZVP_libxc_pbe_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
