set(_li_solv_launcher " ")

if(ENABLE_MPI)
    set(_li_solv_launcher "mpirun -np 2")
endif()

add_integration_test(
  NAME "Li_solvent_Effect"
  LABELS "li_solv"
  COST 100
  LAUNCH_AGENT ${_li_solv_launcher}
  )
