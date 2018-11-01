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
        int sh_mem_size = 10000; //in MB. Virtual memory, does not cost anything if not used
        this->shared_mem = new mrcpp::SharedMemory(mpi::comm_share, sh_mem_size);
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

/** @brief In place addition */
void QMFunction::add(ComplexDouble c, QMFunction inp, double prec) {
    QMFunction tmp(*this);
    this->clear();
    qmfunction::add(*this, 1.0, tmp, c, inp, prec);
    tmp.free();
}

/** @brief In place multiply */
void QMFunction::multiply(QMFunction inp, double prec) {
    QMFunction tmp(*this);
    this->clear();
    qmfunction::multiply(*this, tmp, inp, prec);
    tmp.free();
}

/** @brief In place multiply with scalar */
void QMFunction::rescale(ComplexDouble c) {
    double thrs = mrcpp::MachineZero;
    bool cHasReal = (std::abs(c.real()) > thrs);
    bool cHasImag = (std::abs(c.imag()) > thrs);

    if (cHasReal and cHasImag) {
        QMFunction tmp(*this);
        this->clear();
        qmfunction::add(*this, c, tmp, 0.0, tmp, -1.0);
        tmp.free();
    }
    if (cHasReal and not cHasImag) {
        if (this->hasReal()) this->real().rescale(c.real());
        if (this->hasImag()) this->imag().rescale(c.real());
    }
    if (not cHasReal and not cHasImag) {
        if (this->hasReal()) this->real().setZero();
        if (this->hasImag()) this->imag().setZero();
    }
    if (not cHasReal and cHasImag) {
        double conj = (this->conjugate()) ? -1.0 : 1.0;
        mrcpp::FunctionTree<3> *tmp_re = this->re;
        mrcpp::FunctionTree<3> *tmp_im = this->im;
        if (tmp_re != 0) tmp_re->rescale(c.imag());
        if (tmp_im != 0) tmp_im->rescale(-1.0 * conj * c.imag());
        this->clear();
        this->setReal(tmp_im);
        this->setImag(tmp_re);
    }
}

} //namespace mrchem
