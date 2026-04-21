add_integration_test(
  NAME "cube_parser"
  LABELS "cube_guess;mrchem;cube_parser;short"
  COST 200
  )

set_tests_properties("cube_parser" 
  PROPERTIES 
    ENVIRONMENT "MRCHEM_PYTHON_MODULE=${PROJECT_BINARY_DIR}/lib/${PYMOD_INSTALL_LIBDIR}"
  )

