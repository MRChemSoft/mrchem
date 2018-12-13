#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "PositionOperator.h"
#include "qmfunctions/qmfunction_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

PositionPotential::PositionPotential(int d, const mrcpp::Coord<3> &o)
        : QMPotential(1) {
    auto f = [d, o](const mrcpp::Coord<3> &r) -> double {
        return r[d] - o[d];
    };

    this->func.set(f);
}

void PositionPotential::setup(double prec) {
    if (this->isSetup(prec)) return;
    this->setApplyPrec(prec);

    QMPotential &V = *this;

    if (V.function().hasReal()) MSG_ERROR("Potential not properly cleared");
    if (V.function().hasImag()) MSG_ERROR("Potential not properly cleared");

    Timer timer;
    qmfunction::project(V, this->func, NUMBER::Real, this->apply_prec);
    timer.stop();

    int n = V.function().getNNodes(NUMBER::Total);
    double t = timer.getWallTime();
    Printer::printTree(1, "PositionPotential", n, t);
}

void PositionPotential::clear() {
    this->freeFunctions();  // delete FunctionTree pointers
    this->clearApplyPrec(); // apply_prec = -1
}

} //namespace mrchem
