#include "mrchem.h"

/** @file mathutils.h
 *
 * @brief Collection of stand-alone math related functions
 *
 */

namespace mrchem {
namespace mathutils {

double calc_distance(const double *a, const double *b);

void print_matrix(int level, const DoubleMatrix &M, const std::string &name, int pr = 5);

DoubleMatrix read_matrix_file(const std::string &file);
DoubleMatrix skew_matrix_exp(const DoubleMatrix &A);
ComplexMatrix hermitian_matrix_pow(const ComplexMatrix &A, double b);
ComplexMatrix diagonalize_hermitian_matrix(const ComplexMatrix &A, DoubleVector &diag);
void diagonalize_block(ComplexMatrix &M, ComplexMatrix &U, int nstart, int nsize);

} //namespace mathutils
} //namespace mrchem
