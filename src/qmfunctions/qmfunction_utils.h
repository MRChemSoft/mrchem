#pragma once

#include "mrchem.h"

namespace mrchem {

namespace SPIN { enum type { Paired, Alpha, Beta }; }
namespace NUMBER { enum type { Total, Real, Imag }; }
namespace DENSITY { enum type { Total, Spin, Alpha, Beta }; }


class QMFunction;
typedef std::vector<QMFunction> QMFunctionVector;
class Orbital;
typedef std::vector<Orbital> OrbitalVector;
class Density
typedef std::vector<Density> DensityVector;

namespace qmfunction {
void add(ComplexDouble a, QMFunction &inp_a, bool conj_a,
         ComplexDouble b, QMFunction &inp_b, bool conj_b,
         QMFunction &out, double prec = -1.0);

void linear_combination(ComplexVector &coeff, QMFunctionVector &inp, std::vector<bool> &conj,
                        QMFunction &out,
                        double prec = -1.0);

void linear_combination(ComplexMatrix &U,
                        QMFunctionVector &inp, std::vector<bool> &conj,
                        QMFunctionVector &out, double prec = -1.0);

void pairwise_add(ComplexDouble coeff_a, QMFunctionVector &inp_a, std::vector<bool> &conj_a,
                  ComplexDouble coeff_b, QMFunctionVector &inp_b, std::vector<bool> &conj_b,
                  QMFunctionVector out, double prec);

    void multiply(ComplexDouble coeff, QMFunction &inp, QMFunction &out, bool conj, double prec);

void multiply(double prec, QMFunction &out, double out_coef,
              QMFunction &inp_ab, bool conj_b,
              QMFunction &inp_cd, bool conj_d);

} //namespace qmfunction

} //namespace mrchem
