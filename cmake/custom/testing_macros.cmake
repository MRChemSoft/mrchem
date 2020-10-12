macro(add_integration_test)
  set(oneValueArgs NAME COST LAUNCH_AGENT INITIAL_GUESS)
  set(multiValueArgs LABELS DEPENDS)
  cmake_parse_arguments(_integration_test
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )

  if(_integration_test_LAUNCH_AGENT)
    set(_launch_agent_flag "--launch-agent" "${_integration_test_LAUNCH_AGENT}")
  endif()

  add_test(
    NAME
      ${_integration_test_NAME}
    COMMAND
      ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/test
      --binary $<TARGET_FILE_DIR:mrchem.x>
      --work-dir ${CMAKE_CURRENT_BINARY_DIR}
      ${_launch_agent_flag}
    WORKING_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR}
    )

  set_tests_properties(${_integration_test_NAME}
    PROPERTIES
      LABELS "${_integration_test_LABELS};integration"
    )

  if(_integration_test_COST)
    set_tests_properties(${_integration_test_NAME}
      PROPERTIES
        COST ${_integration_test_COST}
      )
  endif()

  if(_integration_test_DEPENDS)
    set_tests_properties(${_integration_test_NAME}
      PROPERTIES
        DEPENDS ${_integration_test_DEPENDS}
      )
  endif()

  if(_integration_test_INITIAL_GUESS)
    file(
      COPY
        ${_integration_test_INITIAL_GUESS}
      DESTINATION
        ${CMAKE_CURRENT_BINARY_DIR}
      )
  endif()
endmacro()

macro(add_Catch_test)
  set(oneValueArgs NAME COST)
  set(multiValueArgs LABELS DEPENDS REFERENCE_FILES)
  cmake_parse_arguments(_Catch_test
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )

  if (ENABLE_MPI)
      set(_unit_launcher ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1)
  endif()

  add_test(
    NAME
      ${_Catch_test_NAME}
    COMMAND
      ${_unit_launcher}
      $<TARGET_FILE:mrchem-tests>
      [${_Catch_test_NAME}] --success --out
      ${PROJECT_BINARY_DIR}/tests/${_Catch_test_NAME}.log --durations yes
    WORKING_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR}
    )

  set_tests_properties(${_Catch_test_NAME}
    PROPERTIES
      LABELS "${_Catch_test_LABELS};unit"
    )

  if(_Catch_test_COST)
    set_tests_properties(${_Catch_test_NAME}
      PROPERTIES
        COST ${_Catch_test_COST}
      )
  endif()

  if(_Catch_test_DEPENDS)
    set_tests_properties(${_Catch_test_NAME}
      PROPERTIES
        DEPENDS ${_Catch_test_DEPENDS}
      )
  endif()

  if(_Catch_test_REFERENCE_FILES)
    file(
      COPY
        ${_Catch_test_REFERENCE_FILES}
      DESTINATION
        ${CMAKE_CURRENT_BINARY_DIR}
      )
  endif()
endmacro()
