#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "qmfunction_utils.h"
#include "QMFunction.h"

using mrcpp::Timer;
using mrcpp::Printer;
using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

ComplexDouble qmfunction::dot(QMFunction &bra, double bra_conj, QMFunction &ket, double ket_conj) {
    double rr(0.0), ri(0.0), ir(0.0), ii(0.0);
    if (bra.hasReal() and ket.hasReal()) rr = mrcpp::dot(bra.real(), ket.real());
    if (bra.hasReal() and ket.hasImag()) ri = mrcpp::dot(bra.real(), ket.imag());
    if (bra.hasImag() and ket.hasReal()) ir = mrcpp::dot(bra.imag(), ket.real());
    if (bra.hasImag() and ket.hasImag()) ii = mrcpp::dot(bra.imag(), ket.imag());

    double real_part = rr + bra_conj*ket_conj*ii;
    double imag_part = ket_conj*ri - bra_conj*ir;
    return ComplexDouble(real_part, imag_part);
}

void qmfunction::multiply(QMFunction &inp_a, double conj_a,
                          QMFunction &inp_b, double conj_b,
                          QMFunction &out, double prec) {
    { // Real part
        FunctionTreeVector<3> vec;
        if (inp_a.hasReal() and inp_b.hasReal()) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = 1.0;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, inp_a.real());
                mrcpp::build_grid(*tree, inp_b.real());
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.real(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.real());
            }
            vec.push_back(std::make_tuple(1.0, tree));
        }
        if (inp_a.hasImag() and inp_b.hasImag()) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = -1.0*conj_a*conj_b;
            if (prec < 0.0) {
                mrcpp::build_grid(*tree, inp_a.imag());
                mrcpp::build_grid(*tree, inp_b.imag());
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.imag(), 0);
            } else {
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.imag());
            }
            vec.push_back(std::make_tuple(1.0, tree));
        }
        if (vec.size() == 1) {
            out.setReal(&mrcpp::get_func(vec, 0));
            mrcpp::clear(vec, false);
        }
        if (vec.size() == 2) {
            out.alloc(NUMBER::Real);
            mrcpp::build_grid(out.real(), vec);
            mrcpp::add(prec, out.real(), vec, 0);
            mrcpp::clear(vec, true);
        }
    }
    
    { // Imaginary part
        FunctionTreeVector<3> vec;
        if (inp_a.hasReal() and inp_b.hasImag()) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = conj_b;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, inp_a.real());
                mrcpp::build_grid(*tree, inp_b.imag());
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.imag(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.imag());
            }
            vec.push_back(std::make_tuple(1.0, tree));
        }
        if (inp_a.hasImag() and inp_b.hasReal()) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = conj_a;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, inp_a.imag());
                mrcpp::build_grid(*tree, inp_b.real());
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.real(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.real());
            }
            vec.push_back(std::make_tuple(1.0, tree));
        }
        if (vec.size() == 1) {
            out.setImag(&mrcpp::get_func(vec, 0));
            mrcpp::clear(vec, false);
        }
        if (vec.size() == 2) {
            out.alloc(NUMBER::Imag);
            mrcpp::build_grid(out.imag(), vec);
            mrcpp::add(prec, out.imag(), vec, 0);
            mrcpp::clear(vec, true);
        }
    }    
}

void qmfunction::linear_combination(const ComplexVector &c,
                                    QMFunctionVector &inp,
                                    QMFunction &out,
                                    double prec) {
    double thrs = mrcpp::MachineZero;
    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivec;

    for (int i = 0; i < inp.size(); i++) {
        bool cHasReal = (std::abs(c[i].real()) > thrs);
        bool cHasImag = (std::abs(c[i].imag()) > thrs);

        double sign = std::get<0>(inp[i]);
        QMFunction &inp_i = std::get<1>(inp[i]);
        if (cHasReal and inp_i.hasReal()) rvec.push_back(std::make_tuple(      c[i].real(), &inp_i.real()));
        if (cHasImag and inp_i.hasImag()) rvec.push_back(std::make_tuple(-sign*c[i].imag(), &inp_i.imag()));

        if (cHasImag and inp_i.hasReal()) ivec.push_back(std::make_tuple(      c[i].imag(), &inp_i.real()));
        if (cHasReal and inp_i.hasImag()) ivec.push_back(std::make_tuple( sign*c[i].real(), &inp_i.imag()));
    }

    if (rvec.size() > 0) {
        out.alloc(NUMBER::Real);
        if (prec < 0.0) {
            mrcpp::build_grid(out.real(), rvec);
            mrcpp::add(prec, out.real(), rvec, 0);
        } else {
            mrcpp::add(prec, out.real(), rvec);
        }
    }
    if (ivec.size() > 0) {
        out.alloc(NUMBER::Imag);
        if (prec < 0.0) {
            mrcpp::build_grid(out.imag(), ivec);
            mrcpp::add(prec, out.imag(), ivec, 0);
        } else {
            mrcpp::add(prec, out.imag(), ivec);
        }
    }
}

} //namespace mrchem
