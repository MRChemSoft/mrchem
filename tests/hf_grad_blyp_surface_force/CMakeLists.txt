if(ENABLE_MPI)
    set(_hf_grad_lda_launcher_surface "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "HF_molecular_gradient_LDA_surface_integral"
  LABELS "mrchem;gradient;surface_integral"
  COST 200
  LAUNCH_AGENT ${_hf_grad_lda_launcher_surface}
  )
