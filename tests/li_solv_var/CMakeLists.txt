set(_li_solv_var_launcher " ")

if(ENABLE_MPI)
    set(_li_solv_var_launcher "mpirun -np 2")
endif()

add_integration_test(
  NAME "Li_Solv_var"
  LABELS "mrchem;li_solv_var;Li_SOLV_VAR;scf;energy"
  COST 100
  LAUNCH_AGENT ${_li_solv_var_launcher}
  )
