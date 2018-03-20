#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "PowerPotential.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

QMPower::QMPower(double pow, const double *R, double S)
        : QMPotential(1),
          power(pow) {
    Nuclei nucs;
    nucs.push_back("H", R);
    this->func.push_back(nucs[0], S);
}

void QMPower::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);

    if (hasReal()) MSG_ERROR("Potential not properly cleared");
    if (hasImag()) MSG_ERROR("Potential not properly cleared");

    double p = this->power;
    NuclearFunction &nuc_func = this->func;
    auto f = [p, nuc_func] (const double *r) -> double {
        double f_r = nuc_func.evalf(r);
        return std::pow(f_r, p);
    };

    Timer timer;
    alloc(NUMBER::Real);
    mrcpp::build_grid(this->real(), this->func);
    mrcpp::project(this->apply_prec, this->real(), f);
    timer.stop();

    int n = getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(0, "Cubic potential", n, t);
}

void QMPower::clear() {
    free();           // delete FunctionTree pointers
    clearApplyPrec(); // apply_prec = -1
}

} //namespace mrchem
