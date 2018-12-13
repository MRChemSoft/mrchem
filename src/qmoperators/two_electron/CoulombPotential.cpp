#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "CoulombPotential.h"
#include "parallel.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"

using mrcpp::PoissonOperator;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief constructor
 *
 * @param P: MW Poisson operator
 * @param Phi: orbitals defining the operator
 *
 * Density is not spin-separated since this operator requires only total density.
 * This operator will always point to the same OrbitalVector, but the orbitals within
 * the vector can change throughout the calculation. The density and (*this)
 * QMPotential is uninitialized at this point and will be computed at setup.
 */

CoulombPotential::CoulombPotential(PoissonOperator *P)
        : QMPotential(1, mpi::share_coul_pot)
        , density(mpi::share_coul_dens)
        , poisson(P) {}

/** @brief prepare operator for application
 *
 * @param prec: apply precision
 *
 * This will compute the Coulomb potential by application of the Poisson
 * operator to the density. If the density is not available it is computed
 * from the current orbitals (assuming that the orbitals are available).
 * For first-order perturbations the first order density and the Hessian will be
 * computed. In order to make the Hessian available to CoulombOperator, it is stored in the 
 * potential function instead of the zeroth-order potential.
 *
 */
void CoulombPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);
    setupDensity(prec);
    setupPotential(prec);
}

/** @brief clear operator after application
 *
 * This will clear the operator and bring it back to the state after construction.
 * The operator can now be reused after another setup.
 */
void CoulombPotential::clear() {
    QMFunction::freeFunctions();   // delete FunctionTree pointers
    this->density.freeFunctions(); // delete FunctionTree pointers
    clearApplyPrec();              // apply_prec = -1
}

/** @brief compute Coulomb potential
 *
 * @param prec: apply precision
 *
 * This will compute the Coulomb potential by application o the Poisson operator
 * to the precomputed electron density.
 */
void CoulombPotential::setupPotential(double prec) {
    if (this->poisson == nullptr) MSG_ERROR("Poisson operator not initialized");

    PoissonOperator &P = *this->poisson;
    ComplexFunction &V = this->function();
    ComplexFunction &rho = this->density.function();

    if (V.hasReal()) MSG_ERROR("Potential not properly cleared");
    if (V.hasImag()) MSG_ERROR("Potential not properly cleared");

    // Adjust precision by system size
    double abs_prec = prec / this->density.norm();
    bool need_to_apply = not(V.isShared()) or mpi::share_master();

    Timer timer;
    V.alloc(NUMBER::Real);
    if (need_to_apply) mrcpp::apply(abs_prec, V.real(), P, rho.real());
    mpi::share_function(V, 0, 22445, mpi::comm_share);
    timer.stop();

    int n = V.getNNodes(NUMBER::Total);
    double t = timer.getWallTime();
    Printer::printTree(0, "Coulomb potential", n, t);
}

} //namespace mrchem
