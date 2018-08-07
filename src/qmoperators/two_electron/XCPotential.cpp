#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "XCPotential.h"
#include "Orbital.h"
#include "Density.h"

using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Constructor
 *
 * @param[in] F XCFunctional pointer
 * @param[in] Phi Vector of orbitals
 *
 * Based on the order and spin the correct nr. of potential functions is determined.
 * Then the functional is set up for subsequent calculations, fixing some internals of
 * xcfun when F.evalSetup is invoked.
 */
XCPotential::XCPotential(mrdft::XCFunctional *F, OrbitalVector *Phi)
        : QMPotential(1),
          orbitals(Phi),
          functional(F) {
}

/** @brief Prepare the operator for application
 * 
 * @param[in] prec Apply precision
 *
 * Sequence of steps required to compute the XC potentials:
 *
 * 1) Compute density
 * 2) Setup xcfun input functions (gradients etc.)
 * 3) Evaluate xcfun
 * 4) Compute XC potential(s) from xcfun output
 *
 */
void XCPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);
    setupDensity();
    setupPotential(prec);
}

/** @brief Clears all data in the XCPotential object */
void XCPotential::clear() {
    this->energy = 0.0;
    mrcpp::clear(this->potentials, true);
    clearApplyPrec();
}

/** @brief Compute electron density
 *
 * The density is computed on the grid provided by the MRDFT module. The grid
 * is kept as is, e.i. no additional refinement at this point, since the grid
 * size is determined inside the module.
 */
void XCPotential::setupDensity() {
    if (this->functional->hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    OrbitalVector &Phi = *this->orbitals;
    if (this->functional->isSpinSeparated()) {
        Timer time_a;
        FunctionTree<3> &tmp_a = this->functional->getDensity(mrdft::DensityType::Alpha);
        Density rho_a;
        rho_a.setReal(&tmp_a);
        density::compute(-1.0, rho_a, Phi, DENSITY::Alpha);
        time_a.stop();
        Printer::printTree(0, "XC alpha density", rho_a.getNNodes(), time_a.getWallTime());

        Timer time_b;
        FunctionTree<3> &tmp_b = this->functional->getDensity(mrdft::DensityType::Beta);
        Density rho_b;
        rho_b.setReal(&tmp_b);
        density::compute(-1.0, rho_b, Phi, DENSITY::Beta);
        time_b.stop();
        Printer::printTree(0, "XC beta density", rho_b.getNNodes(), time_b.getWallTime());
    } else {
        Timer time_t;
        FunctionTree<3> &tmp_t = this->functional->getDensity(mrdft::DensityType::Total);
        Density rho_t;
        rho_t.setReal(&tmp_t);
        density::compute(-1.0, rho_t, Phi, DENSITY::Total);
        time_t.stop();
        Printer::printTree(0, "XC total density", rho_t.getNNodes(), time_t.getWallTime());
    }
}

/** @brief Compute XC potential(s)
 *
 * @param prec Precision used in refinement of density grid
 *
 * This will invoke a sequence of steps in the XCFunctional to compute the final
 * XC potential(s) that define this operator. Assuming the density has alredy been
 * computed:
 *
 * 1) Setup xcfun input functions (gradients etc.)
 * 2) Evaluate xcfun
 * 3) Compute XC energy by integrating energy density
 * 4) Compute XC potential(s) from xcfun output functions
 * 5) Remove excess grid nodes based on precision
 * 6) Add extra grid nodes based on precision
 * 7) Clear internal functions in XCFunctional (density grid is kept)
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

/** @brief Return FunctionTree for the input density from the XCFunctional
 *
 * @param[in] type Which density to return (alpha, beta or total)
 */
FunctionTree<3> &XCPotential::getDensity(int spin) {
    if (spin == DENSITY::Total) return this->functional->getDensity(mrdft::DensityType::Total);
    if (spin == DENSITY::Alpha) return this->functional->getDensity(mrdft::DensityType::Alpha);
    if (spin == DENSITY::Beta)  return this->functional->getDensity(mrdft::DensityType::Beta);
    MSG_FATAL("Invalid density type");
}

/** @brief Return FunctionTree for the XC spin potential
 *
 * @param[in] type Which spin potential to return (alpha, beta or total)
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
 * @param[in] phi Orbital to which the potential is applied
 *
 * The operator is applied by choosing the correct potential function
 * which is then assigned to the real function part of the operator
 * base-class before the base class function is called.
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
