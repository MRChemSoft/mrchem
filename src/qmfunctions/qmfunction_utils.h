#pragma once

#include "qmfunctions.h"

namespace mrchem {

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
