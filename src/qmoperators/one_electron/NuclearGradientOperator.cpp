#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "NuclearGradientOperator.h"
#include "qmfunctions/qmfunction_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

void NuclearGradientPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);

    QMPotential &V = *this;

    if (V.hasReal()) MSG_ERROR("Potential not properly cleared");
    if (V.hasImag()) MSG_ERROR("Potential not properly cleared");

    Timer timer;
    qmfunction::project(V, this->func, NUMBER::Real, this->apply_prec);
    timer.stop();

    int n = V.getNNodes(NUMBER::Total);
    double t = timer.getWallTime();
    Printer::printTree(0, "Nuclear potential", n, t);
}

void NuclearGradientPotential::clear() {
    free(NUMBER::Total); // delete FunctionTree pointers
    clearApplyPrec();    // apply_prec = -1
}

} //namespace mrchem
