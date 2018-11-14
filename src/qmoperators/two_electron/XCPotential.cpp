#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "XCPotential.h"
#include "parallel.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"

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
        : QMPotential(1, mpi::share_xc_pot)
        , orbitals(Phi)
        , functional(F)
        , energy(0.0) {}

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
    setupDensity(prec);
    setupPotential();
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
void XCPotential::setupDensity(double prec) {
    if (this->functional->hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    OrbitalVector &Phi = *this->orbitals;
    if (this->functional->isSpinSeparated()) {
        Timer time_a;
        Density rho_a(false);
        density::compute(prec, rho_a, Phi, DENSITY::Alpha);
        FunctionTree<3> &func_a = this->functional->getDensity(mrdft::DensityType::Alpha);
        mrcpp::copy_grid(func_a, rho_a.real());
        mrcpp::copy_func(func_a, rho_a.real());
        rho_a.free();
        time_a.stop();
        Printer::printTree(0, "XC alpha density", func_a.getNNodes(), time_a.getWallTime());

        Timer time_b;
        Density rho_b(false);
        density::compute(prec, rho_b, Phi, DENSITY::Beta);
        FunctionTree<3> &func_b = this->functional->getDensity(mrdft::DensityType::Beta);
        mrcpp::copy_grid(func_b, rho_b.real());
        mrcpp::copy_func(func_b, rho_b.real());
        rho_b.free();
        time_b.stop();
        Printer::printTree(0, "XC beta density", func_b.getNNodes(), time_b.getWallTime());

        // Extend to union grid
        while (mrcpp::refine_grid(func_a, func_b)) {}
        while (mrcpp::refine_grid(func_b, func_a)) {}
    } else {
        Timer time_t;
        Density rho_t(false);
        density::compute(prec, rho_t, Phi, DENSITY::Total);
        FunctionTree<3> &func_t = this->functional->getDensity(mrdft::DensityType::Total);
        mrcpp::copy_grid(func_t, rho_t.real());
        mrcpp::copy_func(func_t, rho_t.real());
        rho_t.free();
        time_t.stop();
        Printer::printTree(0, "XC total density", func_t.getNNodes(), time_t.getWallTime());
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
void XCPotential::setupPotential() {
    if (this->functional == nullptr) MSG_ERROR("XCFunctional not initialized");
    if (not this->functional->hasDensity()) MSG_ERROR("XC density not initialized");
    if (this->potentials.size() != 0) MSG_ERROR("Potential not properly cleared");

    this->functional->setup();
    this->functional->evaluate();
    this->energy = this->functional->calcEnergy();
    this->potentials = this->functional->calcPotential();
    this->functional->clear();
}

/** @brief Return FunctionTree for the input density from the XCFunctional
 *
 * @param[in] type Which density to return (alpha, beta or total)
 */
FunctionTree<3> &XCPotential::getDensity(int spin) {
    if (spin == DENSITY::Total) return this->functional->getDensity(mrdft::DensityType::Total);
    if (spin == DENSITY::Alpha) return this->functional->getDensity(mrdft::DensityType::Alpha);
    if (spin == DENSITY::Beta) return this->functional->getDensity(mrdft::DensityType::Beta);
    MSG_FATAL("Invalid density type");
}

/** @brief Return FunctionTree for the XC spin potential
 *
 * @param[in] type Which spin potential to return (alpha, beta or total)
 */
FunctionTree<3> &XCPotential::getPotential(int spin) {
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
