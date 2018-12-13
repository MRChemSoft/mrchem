#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "MomentumOperator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

using mrcpp::DerivativeOperator;
using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

QMMomentum::QMMomentum(int d, DerivativeOperator<3> &D)
        : QMOperator()
        , apply_dir(d)
        , derivative(&D) {}

Orbital QMMomentum::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    if (this->derivative == 0) MSG_ERROR("No derivative operator");

    int dir = this->apply_dir;
    DerivativeOperator<3> &D = *this->derivative;

    Orbital out = inp.paramCopy();

    Timer timer;
    // Calc real part
    if (inp.function().hasImag()) {
        out.function().alloc(NUMBER::Real);
        mrcpp::apply(out.function().real(), D, inp.function().imag(), dir);
        if (inp.conjugate()) out.function().real().rescale(-1.0);
    }
    // Calc imag part
    if (inp.function().hasReal()) {
        out.function().alloc(NUMBER::Imag);
        mrcpp::apply(out.function().imag(), D, inp.function().real(), dir);
        out.function().imag().rescale(-1.0);
    }
    timer.stop();

    int n = out.function().getNNodes(NUMBER::Total);
    double t = timer.getWallTime();
    Printer::printTree(1, "Applied QMMomentum", n, t);

    return out;
}

Orbital QMMomentum::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
