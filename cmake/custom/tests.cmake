#.rst:
#
# autocmake.yml configuration::
#
#   docopt: "--enable-tests Enable tests [default: False]."
#   define: "'-DENABLE_TESTS={0}'.format(arguments['--enable-tests'])"

option(ENABLE_TESTS "Enable test suite" ON)

macro(add_Catch_test _name _labels)
  # _labels is not a list, it's a string... Transform it into a list
  set(labels)
  string(REPLACE ";" " " _labels "${_labels}")
  foreach(_label "${_labels}")
    list(APPEND labels ${_label})
  endforeach()
  unset(_labels)

  add_test(NAME ${_name}
           COMMAND ${PROJECT_BINARY_DIR}/tests/unit-tests [${_name}] --success --out ${PROJECT_BINARY_DIR}/tests/${_name}.log --durations yes
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  if(labels)
    set_tests_properties(${_name} PROPERTIES LABELS "${labels}")
  endif()
endmacro()

if(ENABLE_TESTS)
    enable_testing()
    include(CTest)
    add_subdirectory(tests)
endif()
