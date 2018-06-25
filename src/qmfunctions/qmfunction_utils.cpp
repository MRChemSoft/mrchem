#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"
#include "utils/math_utils.h"
#include "utils/RRMaximizer.h"

#include "qmfunctions.h"

using mrcpp::Timer;
using mrcpp::Printer;
using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

namespace qmfunction {

void add(ComplexDouble a, QMFunction &inp_a, bool conj_a,
         ComplexDouble b, QMFunction &inp_b, bool conj_b,
         QMFunction &out, double prec = -1.0) {
    std::vector<bool> conj = {conj_a, conj_b};
    QMFunctionVector  inp  = {inp_a,  inp_b};
    ComplexVector     coefs(2); coefs(0) = a; coefs(1) = b;   
    linear_combination(coefs, inp, conj, out, prec);
}

void linear_combination(ComplexVector &coeff,
                        QMFunctionVector &inp,
                        std::vector<bool> &conj,
                        QMFunction &out,
                        double prec = -1.0) {

    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivec;
    for (int i = 0; i < inp.size(); i++) {
        bool cHasReal = (std::abs(coeff[i].real()) > thrs);
        bool cHasImag = (std::abs(coeff[i].imag()) > thrs);
        double sign = conj[i] ? : -1.0 : 1.0;

        if (cHasReal and inp[i].hasReal()) rvec.push_back(std::make_tuple(      c[i].real(), &inp[i].real()));
        if (cHasImag and inp[i].hasImag()) rvec.push_back(std::make_tuple(-conj*c[i].imag(), &inp[i].imag()));

        if (cHasImag and inp[i].hasReal()) ivec.push_back(std::make_tuple(      c[i].imag(), &inp[i].real()));
        if (cHasReal and inp[i].hasImag()) ivec.push_back(std::make_tuple( conj*c[i].real(), &inp[i].imag()));
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

void linear_combination(ComplexMatrix &U, QMFunctionVector &inp, std::vector<bool> &conj,
                        QMFunctionVector &out, double prec = -1.0) {

    if (U.rows() != out.size()) MSG_ABORT("Output size mismatch");
    if (U.cols() != inp.size()) MSG_ABORT("Input size mismatch");
    
    for (int i = 0; i < U.rows(); i++) {
        const ComplexVector &c = U.row(i); 
        linear_combination(c, inp, conj, out[i], prec);
    }
}


void pairwise_add(ComplexDouble coeff_a, QMFunctionVector &inp_a, std::vector<bool> &conj_a,
                  ComplexDouble coeff_b, QMFunctionVector &inp_b, std::vector<bool> &conj_b,
                  QMFunctionVector out, double prec) {
    if (inp_a.size() != inp_b.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < inp_a.size(); i++) {
        void add(coeff_a, inp_a[i], conj_a[i], inp_b[i], conj_b[i], out[i], prec);
    }
}

void multiply(ComplexDouble coeff, QMFunction &inp, QMFunction &out, bool conj, double prec) {
    bool cHasReal = (std::abs(coeff.real()) > thrs);
    bool cHasImag = (std::abs(coeff.imag()) > thrs);
    double sign = conj ? : -1.0; 1.0;
    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivac;
    if (cHasReal and inp[i].hasReal()) rvec.push_back(std::make_tuple(        coeff.real(), &inp[i].real()));
    if (cHasImag and inp[i].hasImag()) rvec.push_back(std::make_tuple(-sign * coeff.imag(), &inp[i].imag()));
    if (cHasImag and inp[i].hasReal()) ivec.push_back(std::make_tuple(        coeff.imag(), &inp[i].real()));
    if (cHasReal and inp[i].hasImag()) ivec.push_back(std::make_tuple( sign * coeff.real(), &inp[i].imag()));
    
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

/** @brief out = coef * inp_ab * inp_cd
  *
  * param[in] prec multiplication precision
  * param[in/out] out output function
  * param[in] out_coef multiplicative coefficient 
  * param[in] inp_ab (a(x) + ib(x)) first input function
  * param[in] inp_cd (c(x) + id(x)) second input function
  * param[in] conj_b true if the complex conj of the first function is requested
  * param[in] conj_d true if the complex conj of the second function is requested
  *
  * Both inputs can be interpreted as complex conjugate versions of
  * themselves, by adjusting the respective bool flags.  The
  * multiplicative coefficent is a real number to simplify the overall
  * computation. The multiplication by a complex phase factor can be
  * handled in a separate routine
  *
  */
void multiply(QMFunction &inp_ab, bool conj_b,
              QMFunction &inp_cd, bool conj_d,
              QMFunction &out, double out_coef, double prec) {

    double sign_b = conj_b ? : -1.0 : 1.0;
    double sign_d = conj_d ? : -1.0 : 1.0;
    bool ac = (inp_ab.hasReal() and inp_cd.hasReal());
    bool bd = (inp_ab.hasImag() and inp_cd.hasImag());
    bool ad = (inp_ab.hasReal() and inp_cd.hasImag());
    bool bc = (inp_ab.hasImag() and inp_cd.hasReal());
        
    { // Real part
        FunctionTreeVector<3> vec;
        if (ac) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = out_coef;
            if (prec < 0.0) { // Union grid
                mrcpp::build_grid(*tree, inp_ab.real());
                mrcpp::build_grid(*tree, inp_cd.real());
                mrcpp::multiply(prec, *tree, coef, inp_ab.real(), inp_cd.real(), 0);
            } else { // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_ab.real(), inp_cd.real());
            }
            vec.push_back(std::make_tuple(1.0, tree));
        }
        if (bd) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = -1.0 * out_coef * sign_b * sign_d;
            if (prec < 0.0) {
                mrcpp::build_grid(*tree, inp_ab.imag());
                mrcpp::build_grid(*tree, inp_cd.imag());
                mrcpp::multiply(prec, *tree, coef, inp_ab.imag(), inp_cd.imag(), 0);
            } else {
                mrcpp::multiply(prec, *tree, coef, inp_ab.imag(), inp_cd.imag());
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
        if (ad) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = out_coef * sign_d;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, inp_ab.real());
                mrcpp::build_grid(*tree, inp_cd.imag());
                mrcpp::multiply(prec, *tree, coef, inp_ab.real(), inp_cd.imag(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_ab.real(), inp_cd.imag());
            }
            vec.push_back(std::make_tuple(1.0, tree));
        }
        if (bc) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = out_coef * sign_b;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, inp_ab.imag());
                mrcpp::build_grid(*tree, inp_cd.real());
                mrcpp::multiply(prec, *tree, coef, inp_ab.imag(), inp_cd.real(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_ab.imag(), inp_cd.real());
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

} // namespace qmfunction

} // namespace mrchem

