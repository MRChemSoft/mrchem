if(ENABLE_MPI)
    set(_li_solv_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "Li_solvent_effect"
  LABELS "mrchem;Li;scf;dft;pbe0;unrestricted;pcm;solvent;energy;short"
  COST 100
  LAUNCH_AGENT ${_li_solv_launcher}
  )
