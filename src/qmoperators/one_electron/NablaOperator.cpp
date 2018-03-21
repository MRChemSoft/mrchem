#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "NablaOperator.h"

using mrcpp::DerivativeOperator;
using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

QMNabla::QMNabla(int d, DerivativeOperator<3> &D)
    : QMOperator(),
      apply_dir(d),
      derivative(&D) {
}

Orbital QMNabla::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    if (this->derivative == 0) MSG_ERROR("No derivative operator");

    int dir = this->apply_dir;
    DerivativeOperator<3> &D = *this->derivative;

    Orbital out = inp.paramCopy();

    Timer timer;
    // Calc real part
    if (inp.hasReal()) {
        out.alloc(NUMBER::Real);
        mrcpp::apply(out.real(), D, inp.real(), dir);
    }
    // Calc imag part
    if (inp.hasImag()) {
        out.alloc(NUMBER::Imag);
        mrcpp::apply(out.imag(), D, inp.imag(), dir);
        if (inp.conjugate()) out.real().rescale(-1.0);
    }
    timer.stop();

    int n = out.getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(1, "Applied QMNabla", n, t);

    return out;
}

Orbital QMNabla::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
