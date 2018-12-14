#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "QMPotential.h"
#include "qmfunctions/Orbital.h"

using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/** @brief constructor
 *
 * @param adap: extra refinement in output
 *
 * Initializes the QMFunction with NULL pointers for both real and imaginary part.
 * These must be computed in setup() of derived classes. The initial output grid
 * in application will be a copy of the input orbital but NOT a copy of the
 * potential grid. The argument sets how many extra refinement levels is allowed
 * beyond this initial refinement.
 */
QMPotential::QMPotential(int adap, bool shared)
        : QMFunction(shared)
        , QMOperator()
        , adap_build(adap) {}

/** @brief destructor
 *
 * The potential components should be cleared already at this point using clear().
 */
QMPotential::~QMPotential() {
    if (hasReal()) MSG_ERROR("Potential not cleared");
    if (hasImag()) MSG_ERROR("Potential not cleared");
}

/** @brief apply potential
 *
 * @param inp: orbital on which to apply
 *
 * Computes a new orbital that is the product of this potential and the input
 * orbital. Orbital parameters are copied from input.
 */
Orbital QMPotential::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    Timer timer;
    Orbital out = inp.paramCopy();
    calcRealPart(out, inp, false);
    calcImagPart(out, inp, false);
    timer.stop();

    int n = out.getNNodes(NUMBER::Total);
    double t = timer.getWallTime();
    Printer::printTree(1, "Applied QM potential", n, t);

    return out;
}

/** @brief apply complex cojugate potential
 *
 * @param inp: orbital on which to apply
 *
 * Computes a new orbital that is the product of the complex conjugate of this
 * potential and the input orbital. Orbital parameters are copied from input.
 */
Orbital QMPotential::dagger(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    Timer timer;
    Orbital out = inp.paramCopy();
    calcRealPart(out, inp, true);
    calcImagPart(out, inp, true);
    timer.stop();

    int n = out.getNNodes(NUMBER::Total);
    double t = timer.getWallTime();
    Printer::printTree(1, "Applied QM adjoint potential", n, t);

    return out;
}

/** @brief compute real part of output
 *
 * @param inp: input orbital
 * @param dagger: apply complex conjugate potential
 *
 * Computes the real part of the output orbital. The initial output grid is a
 * copy of the input orbital grid but NOT a copy of the potential grid.
 */
void QMPotential::calcRealPart(Orbital &out, Orbital &inp, bool dagger) {
    int adap = this->adap_build;
    double prec = this->apply_prec;

    if (out.hasReal()) MSG_FATAL("Output not empty");
    if (out.isShared()) MSG_FATAL("Cannot share this function");

    QMFunction &V = *this;
    if (V.hasReal() and inp.hasReal()) {
        double coef = 1.0;
        Orbital tmp = out.paramCopy();
        tmp.alloc(NUMBER::Real);
        mrcpp::copy_grid(tmp.real(), inp.real());
        mrcpp::multiply(prec, tmp.real(), coef, V.real(), inp.real(), adap);
        out.add(1.0, tmp);
    }
    if (V.hasImag() and inp.hasImag()) {
        double coef = -1.0;
        if (dagger) coef *= -1.0;
        if (inp.conjugate()) coef *= -1.0;
        Orbital tmp = out.paramCopy();
        tmp.alloc(NUMBER::Real);
        mrcpp::copy_grid(tmp.real(), inp.imag());
        mrcpp::multiply(prec, tmp.real(), coef, V.imag(), inp.imag(), adap);
        out.add(1.0, tmp);
    }
}

/** @brief compute imaginary part of output
 *
 * @param inp: input orbital
 * @param dagger: apply complex conjugate potential
 *
 * Computes the imaginary part of the output orbital. The initial output grid is a
 * copy of the input orbital grid but NOT a copy of the potential grid.
 */
void QMPotential::calcImagPart(Orbital &out, Orbital &inp, bool dagger) {
    int adap = this->adap_build;
    double prec = this->apply_prec;

    if (out.hasImag()) MSG_FATAL("Output not empty");
    if (out.isShared()) MSG_FATAL("Cannot share this function");

    QMFunction &V = *this;
    if (V.hasReal() and inp.hasImag()) {
        double coef = 1.0;
        if (inp.conjugate()) coef *= -1.0;
        Orbital tmp = out.paramCopy();
        tmp.alloc(NUMBER::Imag);
        mrcpp::copy_grid(tmp.imag(), inp.imag());
        mrcpp::multiply(prec, tmp.imag(), coef, V.real(), inp.imag(), adap);
        out.add(1.0, tmp);
    }
    if (V.hasImag() and inp.hasReal()) {
        double coef = 1.0;
        if (dagger) coef *= -1.0;
        Orbital tmp = out.paramCopy();
        tmp.alloc(NUMBER::Imag);
        mrcpp::copy_grid(tmp.imag(), inp.real());
        mrcpp::multiply(prec, tmp.imag(), coef, V.imag(), inp.real(), adap);
        out.add(1.0, tmp);
    }
}

} //namespace mrchem
