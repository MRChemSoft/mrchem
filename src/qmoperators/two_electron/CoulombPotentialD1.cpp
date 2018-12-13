#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombPotentialD1.h"
#include "qmfunctions/density_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

CoulombPotentialD1::CoulombPotentialD1(mrcpp::PoissonOperator *P, OrbitalVector *Phi)
        : CoulombPotential(P),
          orbitals(Phi) {
}

/** @brief compute electron density
 *
 * @param[in] prec: apply precision
 *
 * This will compute the electron density as the sum of squares of the orbitals.
 */
void CoulombPotentialD1::setupDensity(double prec) {
    if (hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");

    OrbitalVector &Phi = *this->orbitals;
    Density &rho = this->density;

    Timer timer;
    density::compute(prec, rho, Phi, DENSITY::Total);
    timer.stop();
    double t = timer.getWallTime();
    int n = rho.function().getNNodes(NUMBER::Total);
    Printer::printTree(0, "Coulomb density", n, t);
}

} //namespace mrchem
