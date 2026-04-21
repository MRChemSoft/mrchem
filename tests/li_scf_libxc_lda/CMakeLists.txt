if(ENABLE_MPI)
    set(_li_scf_svwn5_libxc_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Li_SCF_Libxc_SVWN5"
  LABELS "mrchem;Li;dft;lda;scf;unrestricted;libxc;energy;short"
  COST 100
  LAUNCH_AGENT ${_li_scf_libxc_svwn5_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
