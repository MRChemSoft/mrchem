if(ENABLE_MPI)
    set(_li_pp_scf_hf_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/psppar.Na ${CMAKE_CURRENT_BINARY_DIR}/psppar.Na COPYONLY)

add_integration_test(
  NAME "Na_PP_SCF_HF"
  LABELS "mrchem;Na;scf;hf;pp;psp;xcfun;medium"
  COST 100
  LAUNCH_AGENT ${_li_pp_scf_hf_launcher}
  )
