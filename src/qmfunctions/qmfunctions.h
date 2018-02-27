#pragma once

#include "mrchem.h"

namespace mrchem {

namespace SPIN { enum type { Paired, Alpha, Beta }; }
namespace NUMBER { enum type { Total, Real, Imag }; }
namespace DENSITY { enum type { Total, Spin, Alpha, Beta }; }


class Orbital;
class OrbitalVector;
namespace orbital {

ComplexDouble dot(Orbital bra, Orbital ket);

Orbital add(ComplexDouble a, Orbital inp_a, ComplexDouble b, Orbital inp_b, double prec = -1.0);
Orbital multiply(Orbital inp_a, Orbital inp_b, double prec = -1.0);

bool compare(const Orbital &orb_a, const Orbital &orb_b);
int compare_occ(const Orbital &orb_a, const Orbital &orb_b);
int compare_spin(const Orbital &orb_a, const Orbital &orb_b);

} //namespace orbital


class Density;
class DensityVector;
namespace density {

} //namespace density


} //namespace mrchem
