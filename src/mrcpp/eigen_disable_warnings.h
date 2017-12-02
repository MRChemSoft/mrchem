#pragma once

// Disable obnoxious warnings from Eigen headers
#if defined __INTEL_COMPILER
//     2196 - routine is both "inline" and "noinline" ("noinline" assumed)
//     ICC 12 generates this warning even without any inline keyword, when defining class methods 'inline' i.e. inside of class body
//     typedef that may be a reference type.
//     279 - controlling expression is constant
//     ICC 12 generates this warning on assert(constant_expression_depending_on_template_params) and frankly this is a legitimate use case.
    #pragma warning push
    #pragma warning disable 2196 //279
//#elif defined __clang__
//     -Wconstant-logical-operand - warning: use of logical && with constant operand; switch to bitwise & or remove constant
//    this is really a stupid warning as it warns on compile-time expressions involving enums
//    #pragma clang diagnostic push
//
//#pragma clang diagnostic ignored "-Wconstant-logical-operand"
#endif

