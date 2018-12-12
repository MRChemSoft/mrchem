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

#include "parallel.h"

#include "QMFunction.h"
#include "qmfunction_utils.h"

using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/** @brief Compute <bra|ket> = int bra^\dag(r) * ket(r) dr.
 *
 *  Notice that the <bra| position is already complex conjugated.
 *
 */
ComplexDouble qmfunction::dot(QMFunction bra, QMFunction ket) {
    ComplexFunction &bra_func = bra.function();
    ComplexFunction &ket_func = ket.function();
    double rr(0.0), ri(0.0), ir(0.0), ii(0.0);
    if (bra_func.hasReal() and ket_func.hasReal()) rr = mrcpp::dot(bra_func.real(), ket_func.real());
    if (bra_func.hasReal() and ket_func.hasImag()) ri = mrcpp::dot(bra_func.real(), ket_func.imag());
    if (bra_func.hasImag() and ket_func.hasReal()) ir = mrcpp::dot(bra_func.imag(), ket_func.real());
    if (bra_func.hasImag() and ket_func.hasImag()) ii = mrcpp::dot(bra_func.imag(), ket_func.imag());

    double bra_conj = (bra.conjugate()) ? -1.0 : 1.0;
    double ket_conj = (ket.conjugate()) ? -1.0 : 1.0;

    double real_part = rr + bra_conj * ket_conj * ii;
    double imag_part = ket_conj * ri - bra_conj * ir;
    return ComplexDouble(real_part, imag_part);
}

/** @brief Deep copy
 *
 * Returns a new function which is a full blueprint copy of the input function.
 * This is achieved by building a new grid for the real and imaginary parts and
 * copying.
 */
void qmfunction::deep_copy(QMFunction &out, QMFunction &inp) {
    ComplexFunction &inp_func = inp.function();
    ComplexFunction &out_func = out.function();

    bool need_to_copy = not(out_func.isShared()) or mpi::share_master();
    if (inp_func.hasReal()) {
        out_func.alloc(NUMBER::Real);
        if (need_to_copy) {
            mrcpp::copy_grid(out_func.real(), inp_func.real());
            mrcpp::copy_func(out_func.real(), inp_func.real());
        }
    }
    if (inp_func.hasImag()) {
        out_func.alloc(NUMBER::Imag);
        if (need_to_copy) {
            mrcpp::copy_grid(out_func.imag(), inp_func.imag());
            mrcpp::copy_func(out_func.imag(), inp_func.imag());
            if (out.conjugate()) out_func.imag().rescale(-1.0);
        }
    }
    mpi::share_function(out_func, 0, 1324, mpi::comm_share);
}

void qmfunction::project(QMFunction &out, std::function<double(const mrcpp::Coord<3> &r)> f, int type, double prec) {
    ComplexFunction &out_func = out.function();
    bool need_to_project = not(out_func.isShared()) or mpi::share_master();
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (not out_func.hasReal()) out_func.alloc(NUMBER::Real);
        if (need_to_project) mrcpp::project<3>(prec, out_func.real(), f);
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (not out_func.hasImag()) out_func.alloc(NUMBER::Imag);
        if (need_to_project) mrcpp::project<3>(prec, out_func.imag(), f);
    }
    mpi::share_function(out_func, 0, 123123, mpi::comm_share);
}

void qmfunction::project(QMFunction &out, mrcpp::RepresentableFunction<3> &f, int type, double prec) {
    ComplexFunction &out_func = out.function();
    bool need_to_project = not(out_func.isShared()) or mpi::share_master();
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (not out_func.hasReal()) out_func.alloc(NUMBER::Real);
        if (need_to_project) mrcpp::build_grid(out_func.real(), f);
        if (need_to_project) mrcpp::project<3>(prec, out_func.real(), f);
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (not out_func.hasImag()) out_func.alloc(NUMBER::Imag);
        if (need_to_project) mrcpp::build_grid(out_func.imag(), f);
        if (need_to_project) mrcpp::project<3>(prec, out_func.imag(), f);
    }
    mpi::share_function(out_func, 0, 132231, mpi::comm_share);
}

/** @brief out = a*inp_a + b*inp_b
 *
 * Recast into linear_combination.
 *
 */
void qmfunction::add(QMFunction &out,
                     ComplexDouble a,
                     QMFunction inp_a,
                     ComplexDouble b,
                     QMFunction inp_b,
                     double prec) {
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
    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivec;

    double thrs = mrcpp::MachineZero;
    for (int i = 0; i < inp.size(); i++) {
        ComplexFunction &func_i = inp[i].function();
        double sign = (inp[i].conjugate()) ? -1.0 : 1.0;

        bool cHasReal = (std::abs(c[i].real()) > thrs);
        bool cHasImag = (std::abs(c[i].imag()) > thrs);

        if (cHasReal and func_i.hasReal()) rvec.push_back(std::make_tuple(c[i].real(), &func_i.real()));
        if (cHasImag and func_i.hasImag()) rvec.push_back(std::make_tuple(-sign * c[i].imag(), &func_i.imag()));

        if (cHasImag and func_i.hasReal()) ivec.push_back(std::make_tuple(c[i].imag(), &func_i.real()));
        if (cHasReal and func_i.hasImag()) ivec.push_back(std::make_tuple(sign * c[i].real(), &func_i.imag()));
    }

    ComplexFunction &out_func = out.function();
    if (rvec.size() > 0 and not out_func.hasReal()) out_func.alloc(NUMBER::Real);
    if (ivec.size() > 0 and not out_func.hasImag()) out_func.alloc(NUMBER::Imag);

    bool need_to_add = not(out_func.isShared()) or mpi::share_master();
    if (need_to_add) {
        if (rvec.size() > 0) {
            if (prec < 0.0) {
                mrcpp::build_grid(out_func.real(), rvec);
                mrcpp::add(prec, out_func.real(), rvec, 0);
            } else {
                mrcpp::add(prec, out_func.real(), rvec);
            }
        } else if (out_func.hasReal()) {
            out_func.real().setZero();
        }
        if (ivec.size() > 0) {
            if (prec < 0.0) {
                mrcpp::build_grid(out_func.imag(), ivec);
                mrcpp::add(prec, out_func.imag(), ivec, 0);
            } else {
                mrcpp::add(prec, out_func.imag(), ivec);
            }
        } else if (out_func.hasImag()) {
            out_func.imag().setZero();
        }
    }
    mpi::share_function(out_func, 0, 9911, mpi::comm_share);
}

/** @brief out = Re(inp_a * inp_b)
 *
 */
void qmfunction::multiply_real(QMFunction &out, QMFunction inp_a, QMFunction inp_b, double prec) {
    double conj_a = (inp_a.conjugate()) ? -1.0 : 1.0;
    double conj_b = (inp_b.conjugate()) ? -1.0 : 1.0;

    ComplexFunction &func_a = inp_a.function();
    ComplexFunction &func_b = inp_b.function();
    ComplexFunction &func_out = out.function();

    bool need_to_multiply = not(func_out.isShared()) or mpi::share_master();

    FunctionTreeVector<3> vec;
    if (func_a.hasReal() and func_b.hasReal()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        if (need_to_multiply) {
            double coef = 1.0;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, func_a.real());
                mrcpp::build_grid(*tree, func_b.real());
                mrcpp::multiply(prec, *tree, coef, func_a.real(), func_b.real(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, func_a.real(), func_b.real());
            }
        }
        vec.push_back(std::make_tuple(1.0, tree));
    }
    if (func_a.hasImag() and func_b.hasImag()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        if (need_to_multiply) {
            double coef = -1.0 * conj_a * conj_b;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, func_a.imag());
                mrcpp::build_grid(*tree, func_b.imag());
                mrcpp::multiply(prec, *tree, coef, func_a.imag(), func_b.imag(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, func_a.imag(), func_b.imag());
            }
        }
        vec.push_back(std::make_tuple(1.0, tree));
    }

    if (vec.size() > 0) {
        if (func_out.hasReal()) {
            if (need_to_multiply) func_out.real().clear();
        } else {
            // All sharing procs must allocate
            func_out.alloc(NUMBER::Real);
        }
    }

    if (need_to_multiply) {
        if (vec.size() == 1) {
            mrcpp::FunctionTree<3> &func_0 = mrcpp::get_func(vec, 0);
            mrcpp::copy_grid(func_out.real(), func_0);
            mrcpp::copy_func(func_out.real(), func_0);
            mrcpp::clear(vec, true);
        } else if (vec.size() == 2) {
            mrcpp::build_grid(func_out.real(), vec);
            mrcpp::add(prec, func_out.real(), vec, 0);
            mrcpp::clear(vec, true);
        } else if (func_out.hasReal()) {
            func_out.real().setZero();
        }
    }
    mpi::share_function(func_out, 0, 9191, mpi::comm_share);
}

/** @brief out = Im(inp_a * inp_b)
 *
 */
void qmfunction::multiply_imag(QMFunction &out, QMFunction inp_a, QMFunction inp_b, double prec) {
    double conj_a = (inp_a.conjugate()) ? -1.0 : 1.0;
    double conj_b = (inp_b.conjugate()) ? -1.0 : 1.0;

    ComplexFunction &func_a = inp_a.function();
    ComplexFunction &func_b = inp_b.function();
    ComplexFunction &func_out = out.function();

    bool need_to_multiply = not(func_out.isShared()) or mpi::share_master();

    FunctionTreeVector<3> vec;
    if (func_a.hasReal() and func_b.hasImag()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        if (need_to_multiply) {
            double coef = conj_b;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, func_a.real());
                mrcpp::build_grid(*tree, func_b.imag());
                mrcpp::multiply(prec, *tree, coef, func_a.real(), func_b.imag(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, func_a.real(), func_b.imag());
            }
        }
        vec.push_back(std::make_tuple(1.0, tree));
    }
    if (func_a.hasImag() and func_b.hasReal()) {
        FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
        if (need_to_multiply) {
            double coef = conj_a;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, func_a.imag());
                mrcpp::build_grid(*tree, func_b.real());
                mrcpp::multiply(prec, *tree, coef, func_a.imag(), func_b.real(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, func_a.imag(), func_b.real());
            }
        }
        vec.push_back(std::make_tuple(1.0, tree));
    }

    if (vec.size() > 0) {
        if (func_out.hasImag()) {
            if (need_to_multiply) func_out.imag().clear();
        } else {
            // All sharing procs must allocate
            func_out.alloc(NUMBER::Imag);
        }
    }

    if (need_to_multiply) {
        if (vec.size() == 1) {
            mrcpp::FunctionTree<3> &func_0 = mrcpp::get_func(vec, 0);
            mrcpp::copy_grid(func_out.imag(), func_0);
            mrcpp::copy_func(func_out.imag(), func_0);
            mrcpp::clear(vec, true);
        } else if (vec.size() == 2) {
            mrcpp::build_grid(func_out.imag(), vec);
            mrcpp::add(prec, func_out.imag(), vec, 0);
            mrcpp::clear(vec, true);
        } else if (func_out.hasImag()) {
            func_out.imag().setZero();
        }
    }
    mpi::share_function(func_out, 0, 9292, mpi::comm_share);
}

} //namespace mrchem
