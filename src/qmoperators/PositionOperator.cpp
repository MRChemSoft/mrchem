#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "PositionOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

QMPosition::QMPosition(int d, const double *o)
        : QMPotential(1) {
    double orig = 0.0;
    if (o != 0) orig = o[d];

    auto f = [d, orig] (const double *r) -> double {
        return r[d] - orig;
    };

    this->func.set(f);
}

void QMPosition::setup(double prec) {
    if (this->isSetup(prec)) return;
    this->setApplyPrec(prec);

    if (this->hasReal()) MSG_ERROR("Potential not properly cleared");
    if (this->hasImag()) MSG_ERROR("Potential not properly cleared");

    Timer timer;
    this->alloc(NUMBER::Real);
    mrcpp::project(this->apply_prec, this->real(), this->func);
    timer.stop();

    int n = this->getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(1, "QMPosition", n, t);
}

void QMPosition::clear() {
    this->free();           // delete FunctionTree pointers
    this->clearApplyPrec(); // apply_prec = -1
}

} //namespace mrchem
