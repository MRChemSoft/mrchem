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

QMFunction::QMFunction(bool share, mrcpp::FunctionTree<3> *r, mrcpp::FunctionTree<3> *i)
        : func_data({0, 0, false})
        , shared_mem(nullptr)
        , re(r)
        , im(i) {
    if (share and mpi::share_size > 1) {
        // Memory size in MB defined in input. Virtual memory, does not cost anything if not used.
        this->shared_mem = new mrcpp::SharedMemory(mpi::comm_share, mpi::shared_memory_size);
    }
}

QMFunction::QMFunction(const QMFunction &func)
        : func_data(func.func_data)
        , shared_mem(nullptr)
        , re(func.re)
        , im(func.im) {
    if (func.isShared()) MSG_FATAL("Cannot shallow copy shared trees");
}

QMFunction &QMFunction::operator=(const QMFunction &func) {
    if (this != &func) {
        if (func.isShared()) MSG_FATAL("Cannot shallow copy shared trees");
        this->func_data = func.func_data;
        this->re = func.re;
        this->im = func.im;
    }
    return *this;
}

QMFunction::~QMFunction() {
    if (this->shared_mem != nullptr) delete this->shared_mem;
}

void QMFunction::alloc(int type) {
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (this->hasReal()) MSG_FATAL("Function not empty");
        this->re = new mrcpp::FunctionTree<3>(*MRA, this->shared_mem);
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (this->hasImag()) MSG_FATAL("Function not empty");
        this->im = new mrcpp::FunctionTree<3>(*MRA, this->shared_mem);
    }
}

void QMFunction::clear(int type) {
    if (type == NUMBER::Real or type == NUMBER::Total) this->re = nullptr;
    if (type == NUMBER::Imag or type == NUMBER::Total) this->im = nullptr;
}

void QMFunction::free(int type) {
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (this->hasReal()) delete this->re;
        this->re = nullptr;
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (this->hasImag()) delete this->im;
        this->im = nullptr;
    }
}

void QMFunction::crop(double prec) {
    if (prec < 0.0) return;
    if (hasReal()) this->real().crop(prec, 1.0, false);
    if (hasImag()) this->imag().crop(prec, 1.0, false);
}

/** @brief Returns the orbital meta data
 *
 * Tree sizes (nChunks) are flushed before return.
 */
FunctionData &QMFunction::getFunctionData() {
    this->func_data.nChunksReal = 0;
    this->func_data.nChunksImag = 0;
    if (this->hasReal()) this->func_data.nChunksReal = real().getNChunksUsed();
    if (this->hasImag()) this->func_data.nChunksImag = imag().getNChunksUsed();
    return this->func_data;
}

int QMFunction::getNNodes(int type) const {
    int nNodes = 0;
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (this->hasReal()) nNodes += this->real().getNNodes();
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (this->hasImag()) nNodes += this->imag().getNNodes();
    }
    return nNodes;
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
    if (this->hasReal()) sq_r = this->real().getSquareNorm();
    if (this->hasImag()) sq_i = this->imag().getSquareNorm();

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
 * MPI: The necessary alloc() routines must be called by
 *      ALL MPI processes if *this function is shared.
 */
void QMFunction::add(ComplexDouble c, QMFunction inp) {
    double thrs = mrcpp::MachineZero;
    bool cHasReal = (std::abs(c.real()) > thrs);
    bool cHasImag = (std::abs(c.imag()) > thrs);

    QMFunction &out = *this;
    if (cHasReal) {
        if (inp.hasReal()) {
            if (not out.hasReal()) {
                if (out.isShared()) MSG_FATAL("Shared function not allocated");
                out.alloc(NUMBER::Real);
                out.real().setZero();
            }
            while (mrcpp::refine_grid(out.real(), inp.real())) {}
            out.real().add(c.real(), inp.real());
        }
        if (inp.hasImag()) {
            if (not out.hasImag()) {
                if (out.isShared()) MSG_FATAL("Shared function not allocated");
                out.alloc(NUMBER::Imag);
                out.imag().setZero();
            }
            double conj = (inp.conjugate()) ? -1.0 : 1.0;
            while (mrcpp::refine_grid(out.imag(), inp.imag())) {}
            out.imag().add(conj * c.real(), inp.imag());
        }
    }
    if (cHasImag) {
        if (inp.hasReal()) {
            if (not out.hasImag()) {
                if (out.isShared()) MSG_FATAL("Shared function not allocated");
                out.alloc(NUMBER::Imag);
                out.imag().setZero();
            }
            while (mrcpp::refine_grid(out.imag(), inp.real())) {}
            out.imag().add(c.imag(), inp.real());
        }
        if (inp.hasImag()) {
            if (not out.hasReal()) {
                if (out.isShared()) MSG_FATAL("Shared function not allocated");
                out.alloc(NUMBER::Real);
                out.real().setZero();
            }
            double conj = (inp.conjugate()) ? -1.0 : 1.0;
            while (mrcpp::refine_grid(out.real(), inp.imag())) {}
            out.real().add(-1.0 * conj * c.imag(), inp.imag());
        }
    }
}

/** @brief In place multiply with scalar.
 *
 * MPI: Cannot be called with a complex argument
 *      if *this function is shared memory.
 */
void QMFunction::rescale(ComplexDouble c) {
    double thrs = mrcpp::MachineZero;
    bool cHasReal = (std::abs(c.real()) > thrs);
    bool cHasImag = (std::abs(c.imag()) > thrs);
    if (cHasReal and cHasImag) {
        if (not this->hasReal()) MSG_FATAL("Real part not allocated");
        if (not this->hasImag()) MSG_FATAL("Imaginary part not allocated");

        // Extend to union grid
        while (mrcpp::refine_grid(this->real(), this->imag())) {}
        while (mrcpp::refine_grid(this->imag(), this->real())) {}

        // Create deep copy
        QMFunction tmp(false);
        tmp.alloc();
        mrcpp::copy_grid(tmp.real(), this->real());
        mrcpp::copy_func(tmp.real(), this->real());
        mrcpp::copy_grid(tmp.imag(), this->imag());
        mrcpp::copy_func(tmp.imag(), this->imag());
        mrcpp::clear_grid(this->real());
        mrcpp::clear_grid(this->imag());
        mrcpp::add(-1.0, this->real(), c.real(), tmp.real(), -c.imag(), tmp.imag());
        mrcpp::add(-1.0, this->imag(), c.real(), tmp.imag(), c.imag(), tmp.real());
        tmp.free();
    }
    if (cHasReal and not cHasImag) {
        if (this->hasReal()) this->real().rescale(c.real());
        if (this->hasImag()) this->imag().rescale(c.real());
    }
    if (not cHasReal and cHasImag) {
        double conj = (this->conjugate()) ? -1.0 : 1.0;
        mrcpp::FunctionTree<3> *tmp_re = this->re;
        mrcpp::FunctionTree<3> *tmp_im = this->im;
        if (tmp_re != nullptr) tmp_re->rescale(c.imag());
        if (tmp_im != nullptr) tmp_im->rescale(-1.0 * conj * c.imag());
        this->setReal(tmp_im);
        this->setImag(tmp_re);
    }
    if (not cHasReal and not cHasImag) {
        if (this->hasReal()) this->real().setZero();
        if (this->hasImag()) this->imag().setZero();
    }
}

} //namespace mrchem
