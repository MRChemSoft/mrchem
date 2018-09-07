macro(add_integration_test)
  set(oneValueArgs NAME COST)
  set(multiValueArgs LABELS DEPENDS REFERENCE_FILES)
  cmake_parse_arguments(_integration_test
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )

  add_test(
    NAME
      ${_integration_test_NAME}
    COMMAND
      ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/test
      --binary=$<TARGET_FILE_DIR:mrchem.x>
      --work-dir=${CMAKE_CURRENT_BINARY_DIR}
      --verbose
      # The input file
    WORKING_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR}
    )

  set_tests_properties(${_integration_test_NAME}
    PROPERTIES
      LABELS "${_integration_test_LABELS}"
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

  if(_integration_test_REFERENCE_FILES)
    file(
      COPY
        ${_integration_test_REFERENCE_FILES}
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

  add_test(
    NAME
      ${_Catch_test_NAME}
    COMMAND
      $<TARGET_FILE:mrchem-tests>
      [${_Catch_test_NAME}] --success --out
      ${PROJECT_BINARY_DIR}/tests/${_Catch_test_NAME}.log --durations yes
    WORKING_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR}
    )

  set_tests_properties(${_Catch_test_NAME}
    PROPERTIES
      LABELS "${_Catch_test_LABELS}"
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
