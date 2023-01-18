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

#include "MRCPP/Parallel"
#include "MRCPP/Printer"

#include "QMFunction.h"
#include "parallel.h"
#include "qmfunctions/qmfunction_utils.h"

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

QMFunction::QMFunction(bool share)
        : func_ptr(std::make_shared<ComplexFunction>(share)) {}

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
    out.conj = not(this->conj);
    return out;
}

void QMFunction::setReal(mrcpp::FunctionTree<3> *tree) {
    if (isShared()) MSG_ABORT("Cannot set in shared function");
    this->func_ptr->re = tree;
}

void QMFunction::setImag(mrcpp::FunctionTree<3> *tree) {
    if (isShared()) MSG_ABORT("Cannot set in shared function");
    this->func_ptr->im = tree;
}

void QMFunction::alloc(int type, mrcpp::MultiResolutionAnalysis<3> *mra) {
    if (mra == nullptr) MSG_ABORT("Invalid argument");
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (hasReal()) MSG_ABORT("Real part already allocated");
        this->func_ptr->re = new mrcpp::FunctionTree<3>(*mra, this->func_ptr->shared_mem_re);
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (hasImag()) MSG_ABORT("Imaginary part already allocated");
        this->func_ptr->im = new mrcpp::FunctionTree<3>(*mra, this->func_ptr->shared_mem_im);
    }
}

void QMFunction::free(int type) {
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (hasReal()) delete this->func_ptr->re;
        this->func_ptr->re = nullptr;
        if (this->func_ptr->shared_mem_re) this->func_ptr->shared_mem_re->clear();
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (hasImag()) delete this->func_ptr->im;
        this->func_ptr->im = nullptr;
        if (this->func_ptr->shared_mem_im) this->func_ptr->shared_mem_im->clear();
    }
}

/** @brief Returns the orbital meta data
 *
 * Tree sizes (nChunks) are flushed before return.
 */
FunctionData &QMFunction::getFunctionData() {
    this->func_ptr->flushFuncData();
    return this->func_ptr->func_data;
}

int QMFunction::getSizeNodes(int type) const {
    int size_mb = 0; // Memory size in kB
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (hasReal()) size_mb += real().getSizeNodes();
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (hasImag()) size_mb += imag().getSizeNodes();
    }
    return size_mb;
}

int QMFunction::getNNodes(int type) const {
    int nNodes = 0;
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (hasReal()) nNodes += real().getNNodes();
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (hasImag()) nNodes += imag().getNNodes();
    }
    return nNodes;
}

int QMFunction::crop(double prec) {
    if (prec < 0.0) return 0;
    bool need_to_crop = not(isShared()) or mpi::share_master();
    int nChunksremoved = 0;
    if (need_to_crop) {
        if (hasReal()) nChunksremoved = real().crop(prec, 1.0, false);
        if (hasImag()) nChunksremoved += imag().crop(prec, 1.0, false);
    }
    mpi::share_function(*this, 0, 7744, mpi::comm_share);
    return nChunksremoved;
}

ComplexDouble QMFunction::integrate() const {
    double int_r = 0.0;
    double int_i = 0.0;
    if (hasReal()) int_r = real().integrate();
    if (hasImag()) int_i = imag().integrate();
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
    double sq_r = -1.0;
    double sq_i = -1.0;
    if (hasReal()) sq_r = real().getSquareNorm();
    if (hasImag()) sq_i = imag().getSquareNorm();

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
    double thrs = mrcpp::MachineZero;
    bool cHasReal = (std::abs(c.real()) > thrs);
    bool cHasImag = (std::abs(c.imag()) > thrs);
    bool outNeedsReal = (cHasReal and inp.hasReal()) or (cHasImag and inp.hasImag());
    bool outNeedsImag = (cHasReal and inp.hasImag()) or (cHasImag and inp.hasReal());

    QMFunction &out = *this;
    bool clearReal(false), clearImag(false);
    if (outNeedsReal and not(out.hasReal())) {
        out.alloc(NUMBER::Real);
        clearReal = true;
    }

    if (outNeedsImag and not(out.hasImag())) {
        out.alloc(NUMBER::Imag);
        clearImag = true;
    }

    bool need_to_add = not(out.isShared()) or mpi::share_master();
    if (need_to_add) {
        if (clearReal) out.real().setZero();
        if (clearImag) out.imag().setZero();
        if (cHasReal and inp.hasReal()) {
            while (mrcpp::refine_grid(out.real(), inp.real())) {}
            out.real().add(c.real(), inp.real());
        }
        if (cHasReal and inp.hasImag()) {
            double conj = (inp.conjugate()) ? -1.0 : 1.0;
            while (mrcpp::refine_grid(out.imag(), inp.imag())) {}
            out.imag().add(conj * c.real(), inp.imag());
        }
        if (cHasImag and inp.hasReal()) {
            while (mrcpp::refine_grid(out.imag(), inp.real())) {}
            out.imag().add(c.imag(), inp.real());
        }
        if (cHasImag and inp.hasImag()) {
            double conj = (inp.conjugate()) ? -1.0 : 1.0;
            while (mrcpp::refine_grid(out.real(), inp.imag())) {}
            out.real().add(-1.0 * conj * c.imag(), inp.imag());
        }
    }
    mpi::share_function(out, 0, 9911, mpi::comm_share);
}

/** @brief In place addition of absolute values.
 *
 * Output is extended to union grid.
 *
 */
void QMFunction::absadd(ComplexDouble c, QMFunction inp) {
    double thrs = mrcpp::MachineZero;
    bool cHasReal = (std::abs(c.real()) > thrs);
    bool cHasImag = (std::abs(c.imag()) > thrs);
    bool outNeedsReal = (cHasReal and inp.hasReal()) or (cHasImag and inp.hasImag());
    bool outNeedsImag = (cHasReal and inp.hasImag()) or (cHasImag and inp.hasReal());

    QMFunction &out = *this;
    bool clearReal(false), clearImag(false);
    if (outNeedsReal and not(out.hasReal())) {
        out.alloc(NUMBER::Real);
        clearReal = true;
    }

    if (outNeedsImag and not(out.hasImag())) {
        out.alloc(NUMBER::Imag);
        clearImag = true;
    }

    bool need_to_add = not(out.isShared()) or mpi::share_master();
    if (need_to_add) {
        if (clearReal) out.real().setZero();
        if (clearImag) out.imag().setZero();
        if (cHasReal and inp.hasReal()) {
            while (mrcpp::refine_grid(out.real(), inp.real())) {}
            out.real().absadd(c.real(), inp.real());
        }
        if (cHasReal and inp.hasImag()) {
            double conj = (inp.conjugate()) ? -1.0 : 1.0;
            while (mrcpp::refine_grid(out.imag(), inp.imag())) {}
            out.imag().absadd(conj * c.real(), inp.imag());
        }
        if (cHasImag and inp.hasReal()) {
            while (mrcpp::refine_grid(out.imag(), inp.real())) {}
            out.imag().absadd(c.imag(), inp.real());
        }
        if (cHasImag and inp.hasImag()) {
            double conj = (inp.conjugate()) ? -1.0 : 1.0;
            while (mrcpp::refine_grid(out.real(), inp.imag())) {}
            out.real().absadd(-1.0 * conj * c.imag(), inp.imag());
        }
    }
    mpi::share_function(out, 0, 9912, mpi::comm_share);
}

/** @brief In place multiply with real scalar. Fully in-place.*/
void QMFunction::rescale(double c) {
    bool need_to_rescale = not(isShared()) or mpi::share_master();
    if (need_to_rescale) {
        if (hasReal()) real().rescale(c);
        if (hasImag()) imag().rescale(c);
    }
    mpi::share_function(*this, 0, 5543, mpi::comm_share);
}

/** @brief In place multiply with complex scalar. Involves a deep copy.*/
void QMFunction::rescale(ComplexDouble c) {
    QMFunction &out = *this;
    QMFunction tmp(isShared());
    qmfunction::deep_copy(tmp, out);
    out.free(NUMBER::Total);
    out.add(c, tmp);
}

} // namespace mrchem
