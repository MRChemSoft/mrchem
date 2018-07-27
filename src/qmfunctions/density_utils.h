#pragma once

#include "mrchem.h"

namespace mrchem {

namespace density {

void compute(double prec, Density &rho, Orbital phi, int spin);
void compute(double prec, Density &rho, OrbitalVector &Phi, int spin);

} //namespace density


} //namespace mrchem
