#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombPotential.h"
#include "Orbital.h"
#include "orbital_utils.h"
#include "density_utils.h"

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
CoulombPotential::CoulombPotential(PoissonOperator *P, OrbitalVector *Phi)
        : QMPotential(1),
          density(), //LUCA: check constructor
          orbitals(Phi),
          poisson(P) {
    density.setReal(new FunctionTree<3>(*MRA));
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
    setupDensity(prec);
    setupPotential(prec);
}

/** @brief clear operator after application
 *
 * This will clear the operator and bring it back to the state after construction.
 * The operator can now be reused after another setup.
 */
void CoulombPotential::clear() {
    QMFunction::free(); // delete FunctionTree pointers
    clearApplyPrec();   // apply_prec = -1
    mrcpp::clear_grid(this->density.real()); // clear MW coefs but keep the grid
}

/** @brief compute electron density
 *
 * @param prec: apply precision
 *
 * This will compute the electron density as the sum of squares of the orbitals.
 */
void CoulombPotential::setupDensity(double prec) {
    if (hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");

    OrbitalVector &Phi = *this->orbitals;
    Density &rho = this->density;

    Timer timer;
    density::compute(-1.0, rho, Phi, DENSITY::Total);
    timer.stop();
    double t = timer.getWallTime();
    int n = rho.real().getNNodes(); //LUCA this should be implemented also for the IMAG part
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
    if (this->poisson == nullptr) MSG_ERROR("Poisson operator not initialized");
    if (hasReal()) MSG_ERROR("Potential not properly cleared");
    if (hasImag()) MSG_ERROR("Potential not properly cleared");

    PoissonOperator &P = *this->poisson;
    QMPotential &V = *this;
    Density &rho = this->density;

    int nPoints = rho.real().getTDim()*rho.real().getKp1_d();
    int inpNodes = rho.getNNodes();

    Timer timer;
    V.alloc(NUMBER::Real);
    mrcpp::apply(prec, V.real(), P, rho.real());
    timer.stop();
    int n = V.getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(0, "Coulomb potential", n, t);

    // Prepare density grid for next iteration
    double abs_prec = prec/rho.real().integrate();
    rho.real().crop(abs_prec, 1.0, false);
    mrcpp::refine_grid(rho.real(), abs_prec);

    int newNodes = rho.real().getNNodes() - inpNodes; //LUCA also imag part here

    println(0, " Coulomb grid size   " << std::setw(21) << inpNodes << std::setw(17) << nPoints*inpNodes);
    println(0, " Coulomb grid change " << std::setw(21) << newNodes << std::setw(17) << nPoints*newNodes);
}

} //namespace mrchem
