#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombPotential.h"

using mrcpp::FunctionTree;
using mrcpp::PoissonOperator;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

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
CoulombPotential::CoulombPotential(PoissonOperator &P, OrbitalVector *Phi)
        : QMPotential(1),
          density(nullptr),
          orbitals(Phi),
          poisson(&P) {
}

FunctionTree<3>& CoulombPotential::getDensity() {
    if (this->density == nullptr) this->density = new Density(*MRA);
    return *this->density;
}

/** @brief prepare operator for application
 *
 * @param prec: apply precision
 *
 * This will compute the Coulomb potential by application of the Poisson
 * operator to the density. If the density is not available it is computed
 * from the current orbitals (assuming that the orbitals are available).
 */
void CoulombPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);

    if (this->density == nullptr) setupDensity(prec);
    setupPotential(prec);
}

/** @brief clear operator after application
 *
 * This will clear the operator and bring it back to the state after construction.
 * The operator can now be reused after another setup.
 */
void CoulombPotential::clear() {
    if (this->density != nullptr) delete this->density;
    QMFunction::free(); // delete FunctionTree pointers
    clearApplyPrec();   // apply_prec = -1
}

/** @brief compute electron density
 *
 * @param prec: apply precision
 *
 * This will compute the electron density as the sum of squares of the orbitals.
 */
void CoulombPotential::setupDensity(double prec) {
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->density == nullptr) this->density = new Density(*MRA);

    OrbitalVector &Phi = *this->orbitals;
    Density &rho = *this->density;

    Timer timer;
    density::calc_density(prec, rho, Phi, SPIN::Paired);
    timer.stop();
    double t = timer.getWallTime();
    int n = rho.getNNodes();
    Printer::printTree(0, "Coulomb density", n, t);
}

/** @brief compute Coulomb potential
 *
 * @param prec: apply precision
 *
 * This will compute the Coulomb potential by application o the Poisson operator
 * to the precomputed electron density.
 */
void CoulombPotential::setupPotential(double prec) {
    if (this->density == 0) MSG_ERROR("Coulomb density not initialized");
    if (this->poisson == 0) MSG_ERROR("Poisson operator not initialized");
    if (hasReal()) MSG_ERROR("Potential not properly cleared");
    if (hasImag()) MSG_ERROR("Potential not properly cleared");

    PoissonOperator &P = *this->poisson;
    QMPotential &V = *this;
    Density &rho = *this->density;

    Timer timer;
    V.alloc(NUMBER::Real);
    mrcpp::apply(prec, V.real(), P, rho);
    timer.stop();
    int n = V.getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(0, "Coulomb potential", n, t);
}

} //namespace mrchem
