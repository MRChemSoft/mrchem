#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "XCPotential.h"
#include "XCPotentialD2.h"
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
XCPotentialD2::XCPotentialD2(mrdft::XCFunctional *F, OrbitalVector *Phi, OrbitalVector *X, OrbitalVector *Y)
        : XCPotential(F, Phi)
        , orbitals_x(X)
        , orbitals_y(Y)
        , pertDensity_t(nullptr)
        , pertDensity_a(nullptr)
        , pertDensity_b(nullptr) {}

XCPotentialD2::~XCPotentialD2() {
    mrcpp::clear(this->potentials, true);
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
void XCPotentialD2::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);
    //    setupDensity(prec);
    setupPerturbedDensity(prec);
    //    setupPotential(prec);
}

/** @brief Clears all data in the XCPotentialD2 object */
void XCPotentialD2::clear() {
    this->energy = 0.0;
    if (this->pertDensity_t != nullptr) this->pertDensity_t;
    if (this->pertDensity_a != nullptr) this->pertDensity_a;
    if (this->pertDensity_b != nullptr) this->pertDensity_b;
    clearApplyPrec();
}

//LUCA This does not work in the case of a non spin separated functional used for an open-shell system!!
void XCPotentialD2::setupPerturbedDensity(double prec) {
    if (this->orbitals_x == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_y == nullptr) MSG_ERROR("Orbitals not initialized");
    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &X = *this->orbitals_x;
    OrbitalVector &Y = *this->orbitals_y;

    if (this->functional->isSpinSeparated()) {
        Density dRho_a = *this->pertDensity_a;
        Density dRho_b = *this->pertDensity_b;
        Timer time_a;
        density::compute(prec, dRho_a, Phi, X, Y, DENSITY::Alpha);
        time_a.stop();
        Printer::printTree(0, "XC perturbed alpha density", dRho_a.function().getNNodes(NUMBER::Total), time_a.getWallTime());

        Timer time_b;
        density::compute(prec, dRho_b, Phi, X, Y, DENSITY::Beta);
        time_b.stop();
        Printer::printTree(0, "XC perturbed beta density", dRho_b.function().getNNodes(NUMBER::Total), time_b.getWallTime());

        // Extend to union grid
        while (mrcpp::refine_grid(dRho_a.function().real(), dRho_b.function().real())) {}
        while (mrcpp::refine_grid(dRho_b.function().real(), dRho_a.function().real())) {}
    } else {
        Density dRho_t = *this->pertDensity_t;
        Timer time_t;
        density::compute(prec, dRho_t, Phi, X, Y, DENSITY::Total);
        time_t.stop();
        Printer::printTree(0, "XC perturbed total density", dRho_t.function().getNNodes(NUMBER::Total), time_t.getWallTime());
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
void XCPotentialD2::setupPotential(double prec) {
    if (this->functional == nullptr) MSG_ERROR("XCFunctional not initialized");
    if (not this->functional->hasDensity()) MSG_ERROR("XC density not initialized");
    if (this->potentials.size() != 0) MSG_ERROR("Potential not properly cleared");

    int inpNodes = this->functional->getNNodes();
    int inpPoints = this->functional->getNPoints();

    this->functional->setup();
    this->functional->evaluate();
    //    this->energy = this->functional->calcEnergy();
    this->potentials = this->functional->calcPotential();
    this->functional->clear();

    int newNodes = this->functional->getNNodes() - inpNodes;
    int newPoints = this->functional->getNPoints() - inpPoints;

    println(0, " XC grid size   " << std::setw(26) << inpNodes << std::setw(17) << inpPoints);
    println(0, " XC grid change " << std::setw(26) << newNodes << std::setw(17) << newPoints);
}

/** @brief Return FunctionTree for the XC spin potential
 *
 * @param[in] type Which spin potential to return (alpha, beta or total)
 */
FunctionTree<3> &XCPotentialD2::getPotential(int orbitalSpin, int densitySpin) {
    int pot_idx = XCPotentialD2::getPotentialIndex(orbitalSpin, densitySpin);
    return mrcpp::get_func(this->potentials, pot_idx);
}

int XCPotentialD2::getPotentialIndex(int orbitalSpin, int densitySpin) {

    int spinFunctional = this->functional->isSpinSeparated() ? 1 : 0;

    //Potential order (spin separated): v_a, v_b, v_aa, v_ab, v_bb
    //Potential order (spin restricted, total density only): v_r, v_rr

    int functional_case = spinFunctional; // 0  1
    functional_case += orbitalSpin << 1;  // 0  2  4  6
    functional_case += densitySpin << 3;  // 0  8 16 24

    //OS Paired Alpha Beta
    //DS Total  Spin  Alpha Beta

    // SF    OS             DS          case   Index
    switch (functional_case) {
        //  0     0 (paired)     0 (total)     0    1
        case (0):
            return 1;
        //  1     0 (paired)     0 (total)     1    not valid
        //  0     1 (alpha )     0 (total)     2    not implemented
        case (2):
            return 1;
        //  1     1 (alpha )     0 (total)     3    not valid
        //  0     2 (beta  )     0 (total)     4    not implemented
        case (4):
            return 1;
        //  1     2 (beta  )     0 (total)     5    not valid
        //  0     3 (unused)     0 (total)     6    not valid
        //  1     3 (unused)     0 (total)     7    not valid
        //  0     0 (paired)     1 (spin )     8    not valid
        //  1     0 (paired)     1 (spin )     9    not valid
        //  0     1 (alpha )     1 (spin )    10    not valid
        //  1     1 (alpha )     1 (spin )    11    not valid
        //  0     2 (beta  )     1 (spin )    12    not valid
        //  1     2 (beta  )     1 (spin )    13    not valid
        //  0     3 (unused)     1 (spin )    14    not valid
        //  1     3 (unused)     1 (spin )    15    not valid
        //  0     0 (paired)     2 (alpha)    16    not implemented
        //  1     0 (paired)     2 (alpha)    17    not implemented
        //  0     1 (alpha )     2 (alpha)    18    not valid
        //  1     1 (alpha )     2 (alpha)    19    2
        case (19):
            return 2;
        //  0     2 (beta  )     2 (alpha)    20    not valid
        //  1     2 (beta  )     2 (alpha)    21    3
        case (21):
            return 3;
        //  0     3 (unused)     2 (alpha)    22    not valid
        //  1     3 (unused)     2 (alpha)    23    not valid
        //  0     0 (paired)     3 (beta )    24    not implemented
        //  1     0 (paired)     3 (beta )    25    not implemented
        //  0     1 (alpha )     3 (beta )    26    not valid
        //  1     1 (alpha )     3 (beta )    27    3
        case (27):
            return 3;
        //  0     2 (beta  )     3 (beta )    28    not valid
        //  1     2 (beta  )     3 (beta )    29    4
        case (29):
            return 4;
        //  0     3 (unused)     3 (beta )    30    not valid
        //  1     3 (unused)     3 (beta )    31    not valid
        default:
            NOT_IMPLEMENTED_ABORT;
    }

    return -1;
}

/** @brief XCPotentialD2 application
 *
 * @param[in] phi Orbital to which the potential is applied
 *
 * The operator is applied by choosing the correct potential function
 * which is then assigned to the real function part of the operator
 * base-class before the base class function is called.
 */
Orbital XCPotentialD2::apply(Orbital phi) {
    bool spinSeparated = this->functional->isSpinSeparated();
    bool totalDens = (pertDensity_t != nullptr);
    bool alphaDens = (pertDensity_a != nullptr);
    bool betaDens = (pertDensity_b != nullptr);

    if (totalDens and (alphaDens or betaDens))
        MSG_ERROR("Total density and spin separated densities are both available");

    FunctionTree<3> *Vrho = new FunctionTree<3>(*MRA);
    mrcpp::FunctionTreeVector<3> components;
    if (not spinSeparated and totalDens) {
        FunctionTree<3> *component = buildComponent(phi.spin(), DENSITY::Total, *pertDensity_t);
        components.push_back(std::make_tuple(1.0, component));
    } else if (spinSeparated) {
        if (totalDens) NOT_IMPLEMENTED_ABORT;
        if (alphaDens) {
            FunctionTree<3> *component = buildComponent(phi.spin(), DENSITY::Alpha, *pertDensity_a);
            components.push_back(std::make_tuple(1.0, component));
        }
        if (betaDens) {
            FunctionTree<3> *component = buildComponent(phi.spin(), DENSITY::Beta, *pertDensity_b);
            components.push_back(std::make_tuple(1.0, component));
        }
    } else {
        NOT_IMPLEMENTED_ABORT;
    }

    mrcpp::build_grid(*Vrho, components); //LUCA just using "add" results in loss of precision.
    mrcpp::add(-1.0, *Vrho, components);
    this->function().setReal(Vrho);
    Orbital Vrhophi = QMPotential::apply(phi);
    this->function().setReal(nullptr);
    delete Vrho;
    mrcpp::clear(components, true);
    return Vrhophi;
}

FunctionTree<3> *XCPotentialD2::buildComponent(int orbital_spin, int density_spin, Density &pert_dens) {
    FunctionTree<3> &dRho = pert_dens.function().real();
    FunctionTree<3> *tmp = new FunctionTree<3>(*MRA);
    FunctionTree<3> &V = getPotential(orbital_spin, density_spin);
    mrcpp::build_grid(*tmp, V);
    mrcpp::build_grid(*tmp, dRho);
    mrcpp::multiply(-1.0, *tmp, 1.0, V, dRho);
    return tmp;
}

} // namespace mrchem
