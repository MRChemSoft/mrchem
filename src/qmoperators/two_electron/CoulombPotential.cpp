#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombPotential.h"

using mrcpp::PoissonOperator;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

CoulombPotential::CoulombPotential(PoissonOperator &P, OrbitalVector &Phi)
        : QMPotential(1),
          density(false, true),
          orbitals(&Phi),
          poisson(&P) {
}

void CoulombPotential::setup(double prec) {
    if (this->isSetup(prec)) return;
    setApplyPrec(prec);

    setupDensity(prec);
    setupPotential(prec);
}

void CoulombPotential::clear() {
    this->free();           // delete FunctionTree pointers
    this->density.free();   // delete FunctionTree pointers
    this->clearApplyPrec(); // apply_prec = -1
}

void CoulombPotential::setupDensity(double prec) {
    if (this->orbitals == 0) MSG_ERROR("Orbitals not initialized");

    Density &rho = this->density;
    OrbitalVector &Phi = *this->orbitals;

    Timer timer;
    density::calc_density(rho, Phi, prec);
    timer.stop();
    double t = timer.getWallTime();
    int n = rho.getNNodes();
    Printer::printTree(0, "Coulomb density", n, t);
}

void CoulombPotential::setupPotential(double prec) {
    if (this->poisson == 0) MSG_ERROR("Poisson operator not initialized");
    if (this->hasReal()) MSG_ERROR("Potential not properly cleared");
    if (this->hasImag()) MSG_ERROR("Potential not properly cleared");

    PoissonOperator &P = *this->poisson;
    QMPotential &V = *this;
    Density &rho = this->density;

    Timer timer;
    V.alloc(NUMBER::Real);
    mrcpp::apply(prec, V.real(), P, rho.total());
    timer.stop();
    int n = V.getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(0, "Coulomb potential", n, t);
}

} //namespace mrchem
