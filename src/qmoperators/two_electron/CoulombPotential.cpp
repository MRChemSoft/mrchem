#include "CoulombPotential.h"
#include "DensityProjector.h"
#include "PoissonOperator.h"
#include "MWConvolution.h"

using namespace std;

void CoulombPotential::calcDensity(Density &rho, OrbitalVector &phi) {
    if (this->orbitals == 0) MSG_ERROR("Orbitals not initialized");

    DensityProjector project(this->apply_prec, this->max_scale);

    Timer timer;
    project(rho, phi);
    timer.stop();
    double t = timer.getWallTime();
    int n = rho.getNNodes();
    TelePrompter::printTree(0, "Coulomb density", n, t);
}

void CoulombPotential::calcPotential(QMPotential &V, Density &rho) {
    if (V.hasReal()) MSG_ERROR("Potential not properly cleared");
    if (V.hasImag()) MSG_ERROR("Potential not properly cleared");

    if (this->poisson == 0) MSG_ERROR("Poisson operator not initialized");
    PoissonOperator &P = *this->poisson;

    MWConvolution<3> apply(this->apply_prec, this->max_scale);

    Timer timer;
    V.allocReal();
    apply(V.real(), P, rho.total());
    timer.stop();
    int n = V.getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(0, "Coulomb potential", n, t);
}
