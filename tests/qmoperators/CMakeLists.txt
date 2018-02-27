add_library(qmoperator-tests OBJECT identity_operator.cpp)

add_Catch_test(identity_operator "identity_operator")
