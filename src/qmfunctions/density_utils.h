#pragma once

#include "mrchem.h"

namespace mrchem {

namespace density {

void compute(double prec, Density &rho, Orbital phi, int spin);
void compute(double prec, Density &rho, OrbitalVector &Phi, int spin);
void project(double prec, Density &rho, mrcpp::GaussExp<3> &dens_exp, int spin);

} //namespace density


} //namespace mrchem
