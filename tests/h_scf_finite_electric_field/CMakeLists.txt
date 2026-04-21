if(ENABLE_MPI)
    set(_h_el_field_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H_SCF_finite_electric_field"
  LABELS "mrchem;h;electric_field;scf;hf;short"
  COST 100
  LAUNCH_AGENT ${_h_el_field_launcher}
  )
