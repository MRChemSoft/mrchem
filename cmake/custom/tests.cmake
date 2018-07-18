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

include(testing_macros)

if(ENABLE_TESTS)
  enable_testing()
  include(CTest)
  add_subdirectory(tests) # This must come last!!
endif()
