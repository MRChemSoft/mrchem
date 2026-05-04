add_integration_test(
  NAME "cube_parser"
  LABELS "guess;mrchem;cube;parser;short"
  COST 200
  )

set_tests_properties("cube_parser" 
  PROPERTIES 
    ENVIRONMENT "MRCHEM_PYTHON_MODULE=${PROJECT_BINARY_DIR}/lib/${PYMOD_INSTALL_LIBDIR}"
  )

