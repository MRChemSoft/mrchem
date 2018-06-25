#pragma once

#include "mrchem.h"

namespace mrchem {

namespace density {

void compute(double prec, Density &rho, Orbital phi, int spin);
void compute(double prec, Density &rho, OrbitalVector &Phi, int spin);
double calc_density_occupancy(int spin_orb, int spin_dens, double occ_orb);
double calc_singlet_density_occupancy(spin_orb_0, spin_orb_1, spin_dens,
                                      occ_orb_0,  occ_orb_1);
double calc_triplet_density_occupancy(spin_orb_0, spin_orb_1, spin_dens,
                                      occ_orb_0, occ_orb_1);
 
//void compute(double prec, Density &rho, OrbitalVector &Phi_0, &Phi_1, int spin);

} //namespace density


} //namespace mrchem
