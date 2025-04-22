/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "MRCPP/Gaussians"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"

#include "mrcpp::ComplexFunction.h"
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
ComplexDouble mrcpp::cplxfunc::dot(mrcpp::ComplexFunction bra, mrcpp::ComplexFunction ket) {
    double rr(0.0), ri(0.0), ir(0.0), ii(0.0);
    if (bra.hasReal() and ket.hasReal()) rr = mrcpp::dot(bra.real(), ket.real());
    if (bra.hasReal() and ket.hasImag()) ri = mrcpp::dot(bra.real(), ket.imag());
    if (bra.hasImag() and ket.hasReal()) ir = mrcpp::dot(bra.imag(), ket.real());
    if (bra.hasImag() and ket.hasImag()) ii = mrcpp::dot(bra.imag(), ket.imag());

    double bra_conj = (bra.conjugate()) ? -1.0 : 1.0;
    double ket_conj = (ket.conjugate()) ? -1.0 : 1.0;

    double real_part = rr + bra_conj * ket_conj * ii;
    double imag_part = ket_conj * ri - bra_conj * ir;
    return ComplexDouble(real_part, imag_part);
}

/** @brief Compute <bra|ket> = int |bra^\dag(r)| * |ket(r)| dr.
 *
 */
ComplexDouble mrcpp::cplxfunc::node_norm_dot(mrcpp::ComplexFunction bra, mrcpp::ComplexFunction ket, bool exact) {
    double rr(0.0), ri(0.0), ir(0.0), ii(0.0);
    if (bra.hasReal() and ket.hasReal()) rr = mrcpp::node_norm_dot(bra.real(), ket.real(), exact);
    if (bra.hasReal() and ket.hasImag()) ri = mrcpp::node_norm_dot(bra.real(), ket.imag(), exact);
    if (bra.hasImag() and ket.hasReal()) ir = mrcpp::node_norm_dot(bra.imag(), ket.real(), exact);
    if (bra.hasImag() and ket.hasImag()) ii = mrcpp::node_norm_dot(bra.imag(), ket.imag(), exact);

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
void mrcpp::cplxfunc::deep_copy(mrcpp::ComplexFunction &out, mrcpp::ComplexFunction &inp) {
    bool need_to_copy = not(out.isShared()) or mrcpp::mpi::share_master();
    if (inp.hasReal()) {
        if (not out.hasReal()) out.alloc(NUMBER::Real);
        if (need_to_copy) {
            mrcpp::copy_grid(out.real(), inp.real());
            mrcpp::copy_func(out.real(), inp.real());
        }
    }
    if (inp.hasImag()) {
        if (not out.hasImag()) out.alloc(NUMBER::Imag);
        if (need_to_copy) {
            mrcpp::copy_grid(out.imag(), inp.imag());
            mrcpp::copy_func(out.imag(), inp.imag());
            if (out.conjugate()) out.imag().rescale(-1.0);
        }
    }
    mrcpp::mpi::share_function(out, 0, 1324, mrcpp::mpi::comm_share);
}

void mrcpp::cplxfunc::project(mrcpp::ComplexFunction &out, std::function<double(const mrcpp::Coord<3> &r)> f, int type, double prec) {
    bool need_to_project = not(out.isShared()) or mrcpp::mpi::share_master();
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (not out.hasReal()) out.alloc(NUMBER::Real);
        if (need_to_project) mrcpp::project<3>(prec, out.real(), f);
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (not out.hasImag()) out.alloc(NUMBER::Imag);
        if (need_to_project) mrcpp::project<3>(prec, out.imag(), f);
    }
    mrcpp::mpi::share_function(out, 0, 123123, mrcpp::mpi::comm_share);
}

void mrcpp::cplxfunc::project(mrcpp::ComplexFunction &out, mrcpp::GaussExp<3> &f, int type, double prec) {
    bool need_to_project = not(out.isShared()) or mrcpp::mpi::share_master();
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (not out.hasReal()) out.alloc(NUMBER::Real);
        if (need_to_project) {
            for (int i = 0; i < f.size(); i++) mrcpp::build_grid(out.real(), *f[i]);
        }
        if (need_to_project) mrcpp::project<3>(prec, out.real(), f);
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (not out.hasImag()) out.alloc(NUMBER::Imag);
        if (need_to_project) {
            for (int i = 0; i < f.size(); i++) mrcpp::build_grid(out.imag(), *f[i]);
        }
        if (need_to_project) mrcpp::project<3>(prec, out.imag(), f);
    }
    mrcpp::mpi::share_function(out, 0, 132231, mrcpp::mpi::comm_share);
}

void mrcpp::cplxfunc::project(mrcpp::ComplexFunction &out, mrcpp::RepresentableFunction<3> &f, int type, double prec) {
    bool need_to_project = not(out.isShared()) or mrcpp::mpi::share_master();
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (not out.hasReal()) out.alloc(NUMBER::Real);
        if (need_to_project) mrcpp::build_grid(out.real(), f);
        if (need_to_project) mrcpp::project<3>(prec, out.real(), f);
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (not out.hasImag()) out.alloc(NUMBER::Imag);
        if (need_to_project) mrcpp::build_grid(out.imag(), f);
        if (need_to_project) mrcpp::project<3>(prec, out.imag(), f);
    }
    mrcpp::mpi::share_function(out, 0, 132231, mrcpp::mpi::comm_share);
}

/** @brief out = a*inp_a + b*inp_b
 *
 * Recast into linear_combination.
 *
 */
void mrcpp::cplxfunc::add(mrcpp::ComplexFunction &out, ComplexDouble a, mrcpp::ComplexFunction inp_a, ComplexDouble b, mrcpp::ComplexFunction inp_b, double prec) {
    ComplexVector coefs(2);
    coefs(0) = a;
    coefs(1) = b;

    mrcpp::ComplexFunctionVector funcs;
    funcs.push_back(inp_a);
    funcs.push_back(inp_b);

    mrcpp::cplxfunc::linear_combination(out, coefs, funcs, prec);
}

/** @brief out = inp_a * inp_b
 *
 */
void mrcpp::cplxfunc::multiply(mrcpp::ComplexFunction &out, mrcpp::ComplexFunction inp_a, mrcpp::ComplexFunction inp_b, double prec, bool absPrec, bool useMaxNorms) {
    multiply_real(out, inp_a, inp_b, prec, absPrec, useMaxNorms);
    multiply_imag(out, inp_a, inp_b, prec, absPrec, useMaxNorms);
}

/** @brief out = c_0*inp_0 + c_1*inp_1 + ... + c_N*inp_N
 *
 */
void mrcpp::cplxfunc::linear_combination(mrcpp::ComplexFunction &out, const ComplexVector &c, ComplexFunctionVector &inp, double prec) {
    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivec;

    double thrs = mrcpp::MachineZero;
    for (int i = 0; i < inp.size(); i++) {
        double sign = (inp[i].conjugate()) ? -1.0 : 1.0;

        bool cHasReal = (std::abs(c[i].real()) > thrs);
        bool cHasImag = (std::abs(c[i].imag()) > thrs);

        if (cHasReal and inp[i].hasReal()) rvec.push_back(std::make_tuple(c[i].real(), &inp[i].real()));
        if (cHasImag and inp[i].hasImag()) rvec.push_back(std::make_tuple(-sign * c[i].imag(), &inp[i].imag()));

        if (cHasImag and inp[i].hasReal()) ivec.push_back(std::make_tuple(c[i].imag(), &inp[i].real()));
        if (cHasReal and inp[i].hasImag()) ivec.push_back(std::make_tuple(sign * c[i].real(), &inp[i].imag()));
    }

    if (rvec.size() > 0 and not out.hasReal()) out.alloc(NUMBER::Real);
    if (ivec.size() > 0 and not out.hasImag()) out.alloc(NUMBER::Imag);

    bool need_to_add = not(out.isShared()) or mrcpp::mpi::share_master();
    if (need_to_add) {
        if (rvec.size() > 0) {
            if (prec < 0.0) {
                mrcpp::build_grid(out.real(), rvec);
                mrcpp::add(prec, out.real(), rvec, 0);
            } else {
                mrcpp::add(prec, out.real(), rvec);
            }
        } else if (out.hasReal()) {
            out.real().setZero();
        }
        if (ivec.size() > 0) {
            if (prec < 0.0) {
                mrcpp::build_grid(out.imag(), ivec);
                mrcpp::add(prec, out.imag(), ivec, 0);
            } else {
                mrcpp::add(prec, out.imag(), ivec);
            }
        } else if (out.hasImag()) {
            out.imag().setZero();
        }
    }
    mrcpp::mpi::share_function(out, 0, 9911, mrcpp::mpi::comm_share);
}

/** @brief out = Re(inp_a * inp_b)
 *
 */
void mrcpp::cplxfunc::multiply_real(mrcpp::ComplexFunction &out, mrcpp::ComplexFunction inp_a, mrcpp::ComplexFunction inp_b, double prec, bool absPrec, bool useMaxNorms) {
    double conj_a = (inp_a.conjugate()) ? -1.0 : 1.0;
    double conj_b = (inp_b.conjugate()) ? -1.0 : 1.0;

    bool need_to_multiply = not(out.isShared()) or mrcpp::mpi::share_master();

    FunctionTreeVector<3> vec;
    if (inp_a.hasReal() and inp_b.hasReal()) {
        auto *tree = new FunctionTree<3>(*MRA);
        if (need_to_multiply) {
            double coef = 1.0;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, inp_a.real());
                mrcpp::build_grid(*tree, inp_b.real());
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.real(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.real(), -1, absPrec, useMaxNorms);
            }
        }
        vec.push_back(std::make_tuple(1.0, tree));
    }
    if (inp_a.hasImag() and inp_b.hasImag()) {
        auto *tree = new FunctionTree<3>(*MRA);
        if (need_to_multiply) {
            double coef = -1.0 * conj_a * conj_b;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, inp_a.imag());
                mrcpp::build_grid(*tree, inp_b.imag());
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.imag(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.imag(), -1, absPrec, useMaxNorms);
            }
        }
        vec.push_back(std::make_tuple(1.0, tree));
    }

    if (vec.size() > 0) {
        if (out.hasReal()) {
            if (need_to_multiply) out.real().clear();
        } else {
            // All sharing procs must allocate
            out.alloc(NUMBER::Real);
        }
    }

    if (need_to_multiply) {
        if (vec.size() == 1) {
            mrcpp::FunctionTree<3> &func_0 = mrcpp::get_func(vec, 0);
            mrcpp::copy_grid(out.real(), func_0);
            mrcpp::copy_func(out.real(), func_0);
            mrcpp::clear(vec, true);
        } else if (vec.size() == 2) {
            mrcpp::build_grid(out.real(), vec);
            mrcpp::add(prec, out.real(), vec, 0);
            mrcpp::clear(vec, true);
        } else if (out.hasReal()) {
            out.real().setZero();
        }
    }
    mrcpp::mpi::share_function(out, 0, 9191, mrcpp::mpi::comm_share);
}

/** @brief out = Im(inp_a * inp_b)
 *
 */
void mrcpp::cplxfunc::multiply_imag(mrcpp::ComplexFunction &out, mrcpp::ComplexFunction inp_a, mrcpp::ComplexFunction inp_b, double prec, bool absPrec, bool useMaxNorms) {
    double conj_a = (inp_a.conjugate()) ? -1.0 : 1.0;
    double conj_b = (inp_b.conjugate()) ? -1.0 : 1.0;

    bool need_to_multiply = not(out.isShared()) or mrcpp::mpi::share_master();

    FunctionTreeVector<3> vec;
    if (inp_a.hasReal() and inp_b.hasImag()) {
        auto *tree = new FunctionTree<3>(*MRA);
        if (need_to_multiply) {
            double coef = conj_b;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, inp_a.real());
                mrcpp::build_grid(*tree, inp_b.imag());
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.imag(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.imag(), -1, absPrec, useMaxNorms);
            }
        }
        vec.push_back(std::make_tuple(1.0, tree));
    }
    if (inp_a.hasImag() and inp_b.hasReal()) {
        auto *tree = new FunctionTree<3>(*MRA);
        if (need_to_multiply) {
            double coef = conj_a;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, inp_a.imag());
                mrcpp::build_grid(*tree, inp_b.real());
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.real(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.real(), -1, absPrec, useMaxNorms);
            }
        }
        vec.push_back(std::make_tuple(1.0, tree));
    }

    if (vec.size() > 0) {
        if (out.hasImag()) {
            if (need_to_multiply) out.imag().clear();
        } else {
            // All sharing procs must allocate
            out.alloc(NUMBER::Imag);
        }
    }

    if (need_to_multiply) {
        if (vec.size() == 1) {
            mrcpp::FunctionTree<3> &func_0 = mrcpp::get_func(vec, 0);
            mrcpp::copy_grid(out.imag(), func_0);
            mrcpp::copy_func(out.imag(), func_0);
            mrcpp::clear(vec, true);
        } else if (vec.size() == 2) {
            mrcpp::build_grid(out.imag(), vec);
            mrcpp::add(prec, out.imag(), vec, 0);
            mrcpp::clear(vec, true);
        } else if (out.hasImag()) {
            out.imag().setZero();
        }
    }
    mrcpp::mpi::share_function(out, 0, 9292, mrcpp::mpi::comm_share);
}

} // namespace mrchem
