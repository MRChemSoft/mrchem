#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "SpinOperator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/qmfunction_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

Orbital QMSpin::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    ComplexDouble coef(0.0, 0.0);
    if (inp.spin() == SPIN::Paired) NOT_IMPLEMENTED_ABORT;
    if (this->D == 0 and inp.spin() == SPIN::Alpha) coef = (0.5, 0.0);
    if (this->D == 0 and inp.spin() == SPIN::Beta) coef = (0.5, 0.0);
    if (this->D == 1 and inp.spin() == SPIN::Alpha) coef = (0.0, 0.5);
    if (this->D == 1 and inp.spin() == SPIN::Beta) coef = (0.0, -0.5);
    if (this->D == 2 and inp.spin() == SPIN::Alpha) coef = (0.5, 0.0);
    if (this->D == 2 and inp.spin() == SPIN::Beta) coef = (-0.5, 0.0);

    Timer timer;
    Orbital out = inp.paramCopy();
    qmfunction::deep_copy(out, inp);
    out.rescale(coef);

    // Flip spin for s_x and s_y
    if (this->D == 0 or this->D == 1) {
        if (inp.spin() == SPIN::Alpha) out.setSpin(SPIN::Beta);
        if (inp.spin() == SPIN::Beta) out.setSpin(SPIN::Alpha);
    }

    timer.stop();
    int n = out.getNNodes(NUMBER::Total);
    double t = timer.getWallTime();
    Printer::printTree(1, "Applied spin operator", n, t);

    return out;
}

Orbital QMSpin::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

} //namespace mrchem
