if(ENABLE_MPI)
    set(_li_scf_pbe0_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Li_SCF_PBE0_Libxc"
  LABELS "mrchem;li_scf_pbe0_libxc;libxc;Li_SCF_PBE0_Libxc;dft;scf;energy"
  COST 100
  LAUNCH_AGENT ${_li_scf_pbe0_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )
