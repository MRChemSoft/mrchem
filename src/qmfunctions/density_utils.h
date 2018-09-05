#pragma once

#include "mrchem.h"

namespace mrchem {

/* The following container classes are defined as derived classes rather than
 * typedefs in order to be able to forward declare them. */

class Orbital;
class OrbitalVector;
class Density;
class DensityChunk final : public std::vector<std::tuple<int, Density> > { };
class DensityVector final : public std::vector<Density> { };

namespace density {

void compute(double prec, Density &rho, Orbital phi, int spin);
void compute(double prec, Density &rho, OrbitalVector &Phi, int spin);
void project(double prec, Density &rho, mrcpp::GaussExp<3> &dens_exp, int spin);

} //namespace density


} //namespace mrchem
