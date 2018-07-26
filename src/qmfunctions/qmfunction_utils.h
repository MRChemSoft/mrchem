#pragma once

#include "mrchem.h"

namespace mrchem {

namespace SPIN { enum type { Paired, Alpha, Beta }; }
namespace NUMBER { enum type { Total, Real, Imag }; }
namespace DENSITY { enum type { Total, Spin, Alpha, Beta }; }


class Orbital;
typedef std::vector<Orbital> OrbitalVector;
typedef std::vector<std::tuple<int, Orbital> > OrbitalChunk;

namespace qmfunction {

} //namespace qmfunction

typedef mrcpp::FunctionTree<3> Density;
typedef mrcpp::FunctionTreeVector<3> DensityVector;

namespace density {

void compute(double prec, Density &rho, Orbital phi, int spin);
void compute(double prec, Density &rho, OrbitalVector &Phi, int spin);

} //namespace density

} //namespace mrchem
