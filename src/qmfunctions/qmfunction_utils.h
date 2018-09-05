#pragma once

#include "mrchem.h"

namespace mrchem {

/* The following container classes are defined as derived classes rather than
 * typedefs in order to be able to forward declare them. */

class QMFunction;
class QMFunctionVector final : public std::vector<std::tuple<double, QMFunction> > { };

namespace qmfunction {
    
ComplexDouble dot(QMFunction &bra, double bra_conj,
                  QMFunction &ket, double ket_conj);
 
void multiply(QMFunction &inp_a, double conj_a,
              QMFunction &inp_b, double conj_b,
              QMFunction &out,   double prec);

void linear_combination(const ComplexVector &c,
                        QMFunctionVector &inp,
                        QMFunction &out,
                        double prec);
 
} //namespace qmfunction

} //namespace mrchem
