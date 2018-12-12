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

#include "MRCPP/Parallel"
#include "MRCPP/Printer"

#include "QMFunction.h"
#include "parallel.h"
#include "qmfunctions/qmfunction_utils.h"

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

QMFunction::QMFunction(bool share)
        : conj(false)
        , func_ptr(std::make_shared<ComplexFunction>(share)) {}

QMFunction::QMFunction(const QMFunction &func)
        : conj(func.conj)
        , func_ptr(func.func_ptr) {}

QMFunction &QMFunction::operator=(const QMFunction &func) {
    if (this != &func) {
        this->conj = func.conj;
        this->func_ptr = func.func_ptr;
    }
    return *this;
}

QMFunction QMFunction::dagger() {
    QMFunction out(*this);
    out.conj = not this->conj;
    return out;
}

void QMFunction::clearFunctions() {
    ComplexFunction &func = this->function();
    bool need_to_clear = not(func.isShared()) or mpi::share_master();
    if (need_to_clear) {
        if (func.hasReal()) func.real().clear();
        if (func.hasImag()) func.imag().clear();
    }
    mpi::share_function(func, 0, 6354, mpi::comm_share);
}

void QMFunction::crop(double prec) {
    if (prec < 0.0) return;
    ComplexFunction &func = this->function();
    bool need_to_crop = not(func.isShared()) or mpi::share_master();
    if (need_to_crop) {
        if (func.hasReal()) func.real().crop(prec, 1.0, false);
        if (func.hasImag()) func.imag().crop(prec, 1.0, false);
    }
    mpi::share_function(func, 0, 7744, mpi::comm_share);
}

ComplexDouble QMFunction::integrate() const {
    const ComplexFunction &func = this->function();
    double int_r = 0.0;
    double int_i = 0.0;
    if (func.hasReal()) int_r = func.real().integrate();
    if (func.hasImag()) int_i = func.imag().integrate();
    return ComplexDouble(int_r, int_i);
}

/** @brief Returns the norm of the orbital */
double QMFunction::norm() const {
    double norm = squaredNorm();
    if (norm > 0.0) norm = std::sqrt(norm);
    return norm;
}

/** @brief Returns the squared norm of the orbital */
double QMFunction::squaredNorm() const {
    const ComplexFunction &func = this->function();
    double sq_r = -1.0;
    double sq_i = -1.0;
    if (func.hasReal()) sq_r = func.real().getSquareNorm();
    if (func.hasImag()) sq_i = func.imag().getSquareNorm();

    double sq_norm = 0.0;
    if (sq_r < 0.0 and sq_i < 0.0) {
        sq_norm = -1.0;
    } else {
        if (sq_r >= 0.0) sq_norm += sq_r;
        if (sq_i >= 0.0) sq_norm += sq_i;
    }
    return sq_norm;
}

/** @brief In place addition.
 *
 * Output is extended to union grid.
 *
 */
void QMFunction::add(ComplexDouble c, QMFunction inp) {
    ComplexFunction &out_func = this->function();
    ComplexFunction &inp_func = inp.function();
    double thrs = mrcpp::MachineZero;
    bool cHasReal = (std::abs(c.real()) > thrs);
    bool cHasImag = (std::abs(c.imag()) > thrs);
    bool outNeedsReal = (cHasReal and inp_func.hasReal()) or (cHasImag and inp_func.hasImag());
    bool outNeedsImag = (cHasReal and inp_func.hasImag()) or (cHasImag and inp_func.hasReal());

    bool clearReal(false), clearImag(false);
    if (outNeedsReal and not(out_func.hasReal())) {
        out_func.alloc(NUMBER::Real);
        clearReal = true;
    }

    if (outNeedsImag and not(out_func.hasImag())) {
        out_func.alloc(NUMBER::Imag);
        clearImag = true;
    }

    bool need_to_add = not(out_func.isShared()) or mpi::share_master();
    if (need_to_add) {
        if (clearReal) out_func.real().setZero();
        if (clearImag) out_func.imag().setZero();
        if (cHasReal and inp_func.hasReal()) {
            while (mrcpp::refine_grid(out_func.real(), inp_func.real())) {}
            out_func.real().add(c.real(), inp_func.real());
        }
        if (cHasReal and inp_func.hasImag()) {
            double conj = (inp.conjugate()) ? -1.0 : 1.0;
            while (mrcpp::refine_grid(out_func.imag(), inp_func.imag())) {}
            out_func.imag().add(conj * c.real(), inp_func.imag());
        }
        if (cHasImag and inp_func.hasReal()) {
            while (mrcpp::refine_grid(out_func.imag(), inp_func.real())) {}
            out_func.imag().add(c.imag(), inp_func.real());
        }
        if (cHasImag and inp_func.hasImag()) {
            double conj = (inp.conjugate()) ? -1.0 : 1.0;
            while (mrcpp::refine_grid(out_func.real(), inp_func.imag())) {}
            out_func.real().add(-1.0 * conj * c.imag(), inp_func.imag());
        }
    }
    mpi::share_function(out_func, 0, 9911, mpi::comm_share);
}

/** @brief In place multiply with real scalar. Fully in-place.*/
void QMFunction::rescale(double c) {
    ComplexFunction &out_func = this->function();
    bool need_to_rescale = not(out_func.isShared()) or mpi::share_master();
    if (need_to_rescale) {
        if (out_func.hasReal()) out_func.real().rescale(c);
        if (out_func.hasImag()) out_func.imag().rescale(c);
    }
    mpi::share_function(out_func, 0, 5543, mpi::comm_share);
}

/** @brief In place multiply with complex scalar. Involves a deep copy.*/
void QMFunction::rescale(ComplexDouble c) {
    QMFunction tmp(this->function().isShared());
    qmfunction::deep_copy(tmp, *this);
    this->freeFunctions();
    this->add(c, tmp);
}

} //namespace mrchem
