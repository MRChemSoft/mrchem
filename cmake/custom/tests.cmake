#.rst:
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--enable-tests Enable integration tests for mrchem.x [default: False]."
#     - "--enable-unit-tests Enable unit tests for libmrchem.a [default: False]."
#   define:
#     - "'-DENABLE_TESTS={0}'.format(arguments['--enable-tests'])"
#     - "'-DENABLE_UNIT_TESTS={0}'.format(arguments['--enable-unit-tests'])"

option(ENABLE_TESTS "Enable test suite" ON)
option(ENABLE_UNIT_TESTS "Enable test suite" ON)

macro(add_Catch_test)
  set(oneValueArgs NAME COST)
  set(multiValueArgs LABELS DEPENDS REFERENCE_FILES)
  cmake_parse_arguments(add_Catch_test
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )

  add_test(
    NAME
      ${add_Catch_test_NAME}
    COMMAND
      $<TARGET_FILE:mrchem-tests>
      [${add_Catch_test_NAME}] --success --out
      ${PROJECT_BINARY_DIR}/tests/${add_Catch_test_NAME}.log --durations yes
    WORKING_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR}
    )

  set_tests_properties(${add_Catch_test_NAME}
    PROPERTIES
      LABELS "${add_Catch_test_LABELS}"
    )

  if(add_Catch_test_COST)
    set_tests_properties(${add_Catch_test_NAME}
      PROPERTIES
        COST ${add_Catch_test_COST}
      )
  endif()

  if(add_Catch_test_DEPENDS)
    set_tests_properties(${add_Catch_test_NAME}
      PROPERTIES
        DEPENDS ${add_Catch_test_DEPENDS}
      )
  endif()

  if(add_Catch_test_REFERENCE_FILES)
    file(
      COPY
        ${add_Catch_test_REFERENCE_FILES}
      DESTINATION
        ${CMAKE_CURRENT_BINARY_DIR}
      )
  endif()
endmacro()
