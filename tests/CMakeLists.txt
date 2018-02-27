add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)

add_executable(unit-tests unit_tests.cpp
    $<TARGET_OBJECTS:qmfunction-tests>
    $<TARGET_OBJECTS:qmoperator-tests>
)
target_link_libraries(unit-tests qmfunctions qmoperators)
