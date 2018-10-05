/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

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

/** @brief Compute <bra|ket> = int bra^\dag(r) * ket(r) dr.
 *
 *  Notice that the <bra| position is already complex conjugated.
 *
 */
ComplexDouble qmfunction::dot(QMFunction bra, QMFunction ket) {
    double rr(0.0), ri(0.0), ir(0.0), ii(0.0);
    if (bra.hasReal() and ket.hasReal()) rr = mrcpp::dot(bra.real(), ket.real());
    if (bra.hasReal() and ket.hasImag()) ri = mrcpp::dot(bra.real(), ket.imag());
    if (bra.hasImag() and ket.hasReal()) ir = mrcpp::dot(bra.imag(), ket.real());
    if (bra.hasImag() and ket.hasImag()) ii = mrcpp::dot(bra.imag(), ket.imag());

    double bra_conj = (bra.conjugate()) ? -1.0 : 1.0;
    double ket_conj = (ket.conjugate()) ? -1.0 : 1.0;

    double real_part = rr + bra_conj*ket_conj*ii;
    double imag_part = ket_conj*ri - bra_conj*ir;
    return ComplexDouble(real_part, imag_part);
}

/** @brief Frees each function in the vector
 *
 * Leaves an empty vector. QMFunctions are freed.
 *
 */
void qmfunction::free(QMFunctionVector &vec) {
    for (int i = 0; i < vec.size(); i++) vec[i].free();
    vec.clear();
}

/** @brief out = a*inp_a + b*inp_b
 *
 * Recast into linear_combination.
 *
 */
void qmfunction::add(QMFunction &out, ComplexDouble a, QMFunction inp_a, ComplexDouble b, QMFunction inp_b, double prec) {
    ComplexVector coefs(2);
    coefs(0) = a;
    coefs(1) = b;

    QMFunctionVector funcs;
    funcs.push_back(inp_a);
    funcs.push_back(inp_b);

    qmfunction::linear_combination(out, coefs, funcs, prec);
}

/** @brief out = inp_a * inp_b
 *
 */
void qmfunction::multiply(QMFunction &out, QMFunction inp_a, QMFunction inp_b, double prec) {
    multiply_real(out, inp_a, inp_b, prec);
    multiply_imag(out, inp_a, inp_b, prec);
}

/** @brief out = c_0*inp_0 + c_1*inp_1 + ... + c_N*inp_N
 *
 */
void qmfunction::linear_combination(QMFunction &out, const ComplexVector &c, QMFunctionVector &inp, double prec) {
    if (out.hasReal()) MSG_ERROR("Output not empty");
    if (out.hasImag()) MSG_ERROR("Output not empty");

    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivec;

    double thrs = mrcpp::MachineZero;
    for (int i = 0; i < inp.size(); i++) {
        double sign = (inp[i].conjugate()) ? -1.0 : 1.0;

        bool cHasReal = (std::abs(c[i].real()) > thrs);
        bool cHasImag = (std::abs(c[i].imag()) > thrs);

        if (cHasReal and inp[i].hasReal()) rvec.push_back(std::make_tuple(      c[i].real(), &inp[i].real()));
        if (cHasImag and inp[i].hasImag()) rvec.push_back(std::make_tuple(-sign*c[i].imag(), &inp[i].imag()));

        if (cHasImag and inp[i].hasReal()) ivec.push_back(std::make_tuple(      c[i].imag(), &inp[i].real()));
        if (cHasReal and inp[i].hasImag()) ivec.push_back(std::make_tuple( sign*c[i].real(), &inp[i].imag()));
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

/** @brief out = Re(inp_a * inp_b)
 *
 */
void qmfunction::multiply_real(QMFunction &out, QMFunction inp_a, QMFunction inp_b, double prec) {
    if (out.hasReal()) MSG_ERROR("Output not empty");

    double conj_a = (inp_a.conjugate()) ? -1.0 : 1.0;
    double conj_b = (inp_b.conjugate()) ? -1.0 : 1.0;

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

/** @brief out = Im(inp_a * inp_b)
 *
 */
void qmfunction::multiply_imag(QMFunction &out, QMFunction inp_a, QMFunction inp_b, double prec) {
    if (out.hasImag()) MSG_ERROR("Output not empty");

    double conj_a = (inp_a.conjugate()) ? -1.0 : 1.0;
    double conj_b = (inp_b.conjugate()) ? -1.0 : 1.0;

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

} //namespace mrchem
