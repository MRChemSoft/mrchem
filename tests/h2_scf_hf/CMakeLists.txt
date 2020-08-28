if(ENABLE_MPI)
    set(_h2_scf_hf_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2_SCF_HartreeFock"
  LABELS "mrchem;h2_scf_hf;H2_SCF_HartreeFock;energy;hartree_fock;scf"
  COST 100
  LAUNCH_AGENT ${_h2_scf_hf_launcher}
  )
