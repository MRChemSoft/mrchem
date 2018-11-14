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
        : QMFunction(shared, nullptr, nullptr)
        , QMOperator()
        , adap_build(adap) {}

/** @brief destructor
 *
 * The potential components should be cleared already at this point using clear().
 */
QMPotential::~QMPotential() {
    if (this->hasReal()) MSG_ERROR("Potential not cleared");
    if (this->hasImag()) MSG_ERROR("Potential not cleared");
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
    FunctionTree<3> *re = calcRealPart(inp, false);
    FunctionTree<3> *im = calcImagPart(inp, false);
    out.set(NUMBER::Real, re);
    out.set(NUMBER::Imag, im);
    timer.stop();

    int n = out.getNNodes();
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
    FunctionTree<3> *re = calcRealPart(inp, true);
    FunctionTree<3> *im = calcImagPart(inp, true);
    out.set(NUMBER::Real, re);
    out.set(NUMBER::Imag, im);
    timer.stop();

    int n = out.getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(1, "Applied QM adjoint potential", n, t);

    return out;
}

/** @brief compute real part of output
 *
 * @param phi: input orbital
 * @param dagger: apply complex conjugate potential
 *
 * Computes the real part of the output orbital. The initial output grid is a
 * copy of the input orbital grid but NOT a copy of the potential grid.
 */
FunctionTree<3> *QMPotential::calcRealPart(Orbital &phi, bool dagger) {
    int adap = this->adap_build;
    double prec = this->apply_prec;

    QMPotential &V = *this;
    FunctionTreeVector<3> vec;

    if (V.hasReal() and phi.hasReal()) {
        double coef = 1.0;
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*tree, phi.real());
        mrcpp::multiply(prec, *tree, coef, V.real(), phi.real(), adap);
        vec.push_back(std::make_tuple(1.0, tree));
    }
    if (V.hasImag() and phi.hasImag()) {
        double coef = -1.0;
        if (dagger) coef *= -1.0;
        if (phi.conjugate()) coef *= -1.0;
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*tree, phi.imag());
        mrcpp::multiply(prec, *tree, coef, V.imag(), phi.imag(), adap);
        vec.push_back(std::make_tuple(1.0, tree));
    }

    FunctionTree<3> *out = 0;
    if (vec.size() == 1) {
        out = &mrcpp::get_func(vec, 0);
        mrcpp::clear(vec, false);
    }
    if (vec.size() == 2) {
        out = new FunctionTree<3>(*MRA);
        mrcpp::build_grid(*out, vec);
        mrcpp::add(-1.0, *out, vec, 0);
        mrcpp::clear(vec, true);
    }
    return out;
}

/** @brief compute imaginary part of output
 *
 * @param phi: input orbital
 * @param dagger: apply complex conjugate potential
 *
 * Computes the imaginary part of the output orbital. The initial output grid is a
 * copy of the input orbital grid but NOT a copy of the potential grid.
 */
FunctionTree<3> *QMPotential::calcImagPart(Orbital &phi, bool dagger) {
    int adap = this->adap_build;
    double prec = this->apply_prec;

    QMPotential &V = *this;
    FunctionTreeVector<3> vec;

    if (V.hasReal() and phi.hasImag()) {
        double coef = 1.0;
        if (phi.conjugate()) coef *= -1.0;
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*tree, phi.imag());
        mrcpp::multiply(prec, *tree, coef, V.real(), phi.imag(), adap);
        vec.push_back(std::make_tuple(1.0, tree));
    }
    if (V.hasImag() and phi.hasReal()) {
        double coef = 1.0;
        if (dagger) coef *= -1.0;
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*tree, phi.real());
        mrcpp::multiply(prec, *tree, coef, V.imag(), phi.real(), adap);
        vec.push_back(std::make_tuple(1.0, tree));
    }

    FunctionTree<3> *out = 0;
    if (vec.size() == 1) {
        out = &mrcpp::get_func(vec, 0);
        mrcpp::clear(vec, false);
    }
    if (vec.size() == 2) {
        out = new FunctionTree<3>(*MRA);
        mrcpp::build_grid(*out, vec);
        mrcpp::add(-1.0, *out, vec, 0);
        mrcpp::clear(vec, true);
    }
    return out;
}

} //namespace mrchem
