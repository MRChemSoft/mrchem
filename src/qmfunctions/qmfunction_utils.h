#pragma once

#include "mrchem.h"

namespace mrchem {

namespace SPIN { enum type { Paired, Alpha, Beta }; }
namespace NUMBER { enum type { Total, Real, Imag }; }
namespace DENSITY { enum type { Total, Spin, Alpha, Beta }; }


class QMFunction;
typedef std::vector<std::tuple<double,QMFunction> > QMFunctionVector;

class Orbital;
typedef std::vector<Orbital> OrbitalVector;
typedef std::vector<std::tuple<int, Orbital> > OrbitalChunk;

class Density;
typedef std::vector<Density> DensityVector;
typedef std::vector<std::tuple<int, Density> > DensityChunk;

namespace qmfunction {
    
ComplexDouble dot(QMFunction &bra, double bra_conj, QMFunction &ket, double ket_conj);
 
void multiply(QMFunction &inp_a, double conj_a,
              QMFunction &inp_b, double conj_b,
              QMFunction &out,   double prec);

void linear_combination(const ComplexVector &c,
                        QMFunctionVector &inp,
                        QMFunction &out,
                        double prec);
 
} //namespace qmfunction

} //namespace mrchem
