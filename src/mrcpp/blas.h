#pragma once

#ifdef HAVE_BLAS
#include "FCMangle.h"

// Fortran function template
extern "C" {
void FC_dgemm(char * transa, char * transb, int * m, int * n, int * k, double * alpha, const double *a, int * lda, const double *b, int * ldb, double * beta, double *c, int * ldc);
}

// C wrapper
inline void dgemm(char transa, char transb, int m, int n, int k, double alpha, const double *a, int lda, const double *b, int ldb, double beta, double *c, int ldc) {
  FC_dgemm(&transa, &transb, &m, &n, &k, &alpha, (double *) a, &lda, (double *) b, &ldb, &beta, c, &ldc);
}

// Fortran function template
extern "C" {
double FC_ddot(int * n, double * dx, int * incx, double * dy, int * incy);
}

// C wrapper
inline double ddot(int n, const double *dx, int incx, const double *dy, int incy) {
  // Fortran takes values by reference
  return FC_ddot(&n, (double *) dx, &incx, (double *) dy, &incy);
}
#endif
