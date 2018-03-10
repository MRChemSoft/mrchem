#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "QMPotential.h"
#include "Orbital.h"

using mrcpp::FunctionTreeVector;
using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

QMPotential::QMPotential(int adap)
        : QMFunction(0, 0),
          QMOperator(),
          adap_build(adap) {
}

QMPotential::~QMPotential() {
    if (this->hasReal()) MSG_ERROR("Potential not cleared");
    if (this->hasImag()) MSG_ERROR("Potential not cleared");
}

Orbital QMPotential::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    Timer timer;
    Orbital out = inp.paramCopy();
    FunctionTree<3> *re = calcRealPart(inp, false);
    FunctionTree<3> *im = calcImagPart(inp, false);
    out.setReal(re);
    out.setImag(im);
    timer.stop();

    int n = out.getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(1, "Applied QM potential", n, t);

    return out;
}

Orbital QMPotential::dagger(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    Timer timer;
    Orbital out = inp.paramCopy();
    FunctionTree<3> *re = calcRealPart(inp, true);
    FunctionTree<3> *im = calcImagPart(inp, true);
    out.setReal(re);
    out.setImag(im);
    timer.stop();

    int n = out.getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(1, "Applied QM adjoint potential", n, t);

    return out;
}

FunctionTree<3>* QMPotential::calcRealPart(Orbital &phi, bool dagger) {
    int adap = this->adap_build;
    double prec = this->apply_prec;

    QMPotential &V = *this;
    FunctionTreeVector<3> vec;

    if (V.hasReal() and phi.hasReal()) {
        double coef = 1.0;
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*tree, phi.real());
        mrcpp::multiply(prec, *tree, coef, V.real(), phi.real(), adap);
        vec.push_back(tree);
    }
    if (V.hasImag() and phi.hasImag()) {
        double coef = -1.0;
        if (dagger) coef *= -1.0;
        if (phi.conjugate()) coef *= -1.0;
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*tree, phi.imag());
        mrcpp::multiply(prec, *tree, coef, V.imag(), phi.imag(), adap);
        vec.push_back(tree);
    }

    FunctionTree<3> *out = 0;
    if (vec.size() == 1) {
        out = vec[0];
        vec.clear(false);
    }
    if (vec.size() == 2) {
        out = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*out, vec);
        mrcpp::add(-1.0, *out, vec, 0);
        vec.clear(true);
    }
    return out;
}

FunctionTree<3>* QMPotential::calcImagPart(Orbital &phi, bool dagger) {
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
        vec.push_back(tree);
    }
    if (V.hasImag() and phi.hasReal()) {
        double coef = 1.0;
        if (dagger) coef *= -1.0;
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*tree, phi.real());
        mrcpp::multiply(prec, *tree, coef, V.imag(), phi.real(), adap);
        vec.push_back(tree);
    }

    FunctionTree<3> *out = 0;
    if (vec.size() == 1) {
        out = vec[0];
        vec.clear(false);
    }
    if (vec.size() == 2) {
        out = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*out, vec);
        mrcpp::add(-1.0, *out, vec, 0);
        vec.clear(true);
    }
    return out;
}

} //namespace mrchem
