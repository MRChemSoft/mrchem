if(ENABLE_MPI)
    set(_hf_grad_lda_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "HF_molecular_gradient_LDA"
  LABELS "hf_grad_lda;mrchem;gradient;HF_molecular_gradient_LDA"
  COST 200
  LAUNCH_AGENT ${_hf_grad_lda_launcher}
  )
