#include "mrchem.h"

namespace mrchem {
namespace mathutils {

double calc_distance(const double *a, const double *b);

DoubleMatrix read_matrix_file(const std::string &file);
DoubleMatrix skew_matrix_exp(const DoubleMatrix &A);
ComplexMatrix hermitian_matrix_pow(const ComplexMatrix &A, double b);
ComplexMatrix diagonalize_hermitian_matrix(const ComplexMatrix &A, DoubleVector &diag);
void diagonalize_block(ComplexMatrix &M, ComplexMatrix &U, int nstart, int nsize);

} //namespace mathutils
} //namespace mrchem
