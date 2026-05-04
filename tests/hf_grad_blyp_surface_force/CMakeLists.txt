if(ENABLE_MPI)
    set(_hf_grad_blyp_surface_force "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "HF_gradient_BLYP_surface_force"
  LABELS "mrchem;HF;dft;gga;blyp;gradient;surface;force;xcfun;short"
  COST 200
  LAUNCH_AGENT ${_hf_grad_blyp_surface_force}
  )
