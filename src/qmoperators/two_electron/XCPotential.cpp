#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "XCPotential.h"
#include "Orbital.h"

using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief constructor
 *
 * @param[in] k order of the operator
 * @param[in] F XCFunctional pointer
 * @param[in] Phi vector of orbitals
 *
 * Based on the order and spin the correct nr. of potential functions is determined
 * Then the functional is set up for subsequent calculations, fixing some internals of
 * xcfun when F.evalSetup is invoked.
 *
 */
XCPotential::XCPotential(mrdft::XCFunctional *F, OrbitalVector *Phi)
        : QMPotential(1),
          orbitals(Phi),
          functional(F) {
}

/** @brief setup the XCPotential
 * 
 * @param[in] prec apply precision
 *
 * Sequence of steps required to compute the XC potentials The most
 * important steps are evaluateXCFunctional and calcPotential where
 * the functional derivatives are computed and the potential assembled
 * respectively
 *
 */
void XCPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);
    setupDensity();
    setupPotential(prec);
}

/** @brief clears all data in the XCPotential object
 *
 */
void XCPotential::clear() {
    this->energy = 0.0;
    mrcpp::clear(this->potentials, true);
    clearApplyPrec();
}

void XCPotential::setupDensity() {
    if (this->functional->hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");

    OrbitalVector &Phi = *this->orbitals;
    if (this->functional->isSpinSeparated()) {
        Timer time_a;
        Density &rho_a = this->functional->getDensity(mrdft::DensityType::Alpha);
        density::compute(-1.0, rho_a, Phi, DENSITY::Alpha);
        time_a.stop();
        Printer::printTree(0, "XC alpha density", rho_a.getNNodes(), time_a.getWallTime());

        Timer time_b;
        Density &rho_b = this->functional->getDensity(mrdft::DensityType::Beta);
        density::compute(-1.0, rho_b, Phi, DENSITY::Beta);
        time_b.stop();
        Printer::printTree(0, "XC beta density", rho_b.getNNodes(), time_b.getWallTime());
    } else {
        Timer time_t;
        Density &rho_t = this->functional->getDensity(mrdft::DensityType::Total);
        density::compute(-1.0, rho_t, Phi, DENSITY::Total);
        time_t.stop();
        Printer::printTree(0, "XC total density", rho_t.getNNodes(), time_t.getWallTime());
    }
}

/** @brief compute XC potential
 *
 * @param prec: apply precision (not used atm)
 *
 */
void XCPotential::setupPotential(double prec) {
    if (this->functional == nullptr) MSG_ERROR("XCFunctional not initialized");
    if (not this->functional->hasDensity()) MSG_ERROR("XC density not initialized");
    if (this->potentials.size() != 0) MSG_ERROR("Potential not properly cleared");

    int inpNodes = this->functional->getNNodes();
    int inpPoints = this->functional->getNPoints();

    this->functional->setup();
    this->functional->evaluate();
    this->energy = this->functional->calcEnergy();
    this->potentials = this->functional->calcPotential();
    this->functional->pruneGrid(prec);
    this->functional->refineGrid(prec);
    this->functional->clear();

    int newNodes = this->functional->getNNodes() - inpNodes;
    int newPoints = this->functional->getNPoints() - inpPoints;

    println(0, " XC grid size   " << std::setw(26) << inpNodes << std::setw(17) << inpPoints);
    println(0, " XC grid change " << std::setw(26) << newNodes << std::setw(17) << newPoints);
}

Density &XCPotential::getDensity(int spin) {
    if (spin == DENSITY::Total) return this->functional->getDensity(mrdft::DensityType::Total);
    if (spin == DENSITY::Alpha) return this->functional->getDensity(mrdft::DensityType::Alpha);
    if (spin == DENSITY::Beta)  return this->functional->getDensity(mrdft::DensityType::Beta);
    MSG_FATAL("Invalid density type");
}

/** @brief fetches the correct index for the potential function to use
 *
 * @param[in] spin of orbital on which to apply
 *
 * Based on the orbital spin, and whether the functional is spin
 * separated, the correct potential index is selected.
 *
 */
FunctionTree<3>& XCPotential::getPotential(int spin) {
    bool spinFunctional = this->functional->isSpinSeparated();
    int pot_idx = -1;
    if (spinFunctional and spin == SPIN::Alpha) {
        pot_idx = 0;
    } else if (spinFunctional and spin == SPIN::Beta) {
        pot_idx = 1;
    } else if (not spinFunctional) {
        pot_idx = 0;
    } else {
        NOT_IMPLEMENTED_ABORT;
    }
    return mrcpp::get_func(this->potentials, pot_idx);
}

/** @brief XCPotential application
 *
 * The operator is applied by choosing the correct potential function
 * which is then assigned to the real function part of the operator
 * base-class before the base class function is called.
 *
 * @param[in] phi orbital to which the potential is applied.
 *
 */
Orbital XCPotential::apply(Orbital phi) {
    if (this->hasImag()) MSG_ERROR("Imaginary part of XC potential non-zero");

    FunctionTree<3> &V = getPotential(phi.spin());
    this->setReal(&V);
    Orbital Vphi = QMPotential::apply(phi); 
    this->setReal(0);

    return Vphi;
}

} //namespace mrchem
