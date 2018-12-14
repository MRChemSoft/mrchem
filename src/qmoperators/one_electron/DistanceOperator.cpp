#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "DistanceOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

DistancePotential::DistancePotential(double pow, const mrcpp::Coord<3> &R, double S)
        : QMPotential(1)
        , power(pow) {
    Nuclei nucs;
    nucs.push_back("H", R);
    this->func.push_back(nucs[0], S);
}

void DistancePotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);

    QMPotential &V = *this;

    if (V.hasReal()) MSG_ERROR("Potential not properly cleared");
    if (V.hasImag()) MSG_ERROR("Potential not properly cleared");
    if (V.isShared()) MSG_FATAL("Cannot share this operator");

    double p = this->power;
    NuclearFunction &nuc_func = this->func;
    auto f = [p, nuc_func](const mrcpp::Coord<3> &r) -> double {
        double f_r = nuc_func.evalf(r);
        return std::pow(f_r, p);
    };

    Timer timer;
    V.alloc(NUMBER::Real);
    mrcpp::build_grid(V.real(), this->func);
    mrcpp::project<3>(this->apply_prec, V.real(), f);
    timer.stop();

    int n = V.getNNodes(NUMBER::Total);
    double t = timer.getWallTime();
    Printer::printTree(0, "Distance potential", n, t);
}

void DistancePotential::clear() {
    free(NUMBER::Total); // delete FunctionTree pointers
    clearApplyPrec();    // apply_prec = -1
}

} //namespace mrchem
