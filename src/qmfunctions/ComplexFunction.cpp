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

#include "ComplexFunction.h"

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

ComplexFunction::ComplexFunction(bool share)
        : func_data({0, 0, share})
        , shared_mem(nullptr)
        , re(nullptr)
        , im(nullptr) {
    if (this->isShared() and mpi::share_size > 1) {
        // Memory size in MB defined in input. Virtual memory, does not cost anything if not used.
        this->shared_mem = new mrcpp::SharedMemory(mpi::comm_share, mpi::shared_memory_size);
    }
}
ComplexFunction::~ComplexFunction() {
    if (this->shared_mem != nullptr) delete this->shared_mem;
    if (this->re != nullptr) delete this->re;
    if (this->im != nullptr) delete this->im;
}

void ComplexFunction::setReal(mrcpp::FunctionTree<3> *tree) {
    if (isShared()) MSG_FATAL("Cannot set in shared function");
    this->re = tree;
}

void ComplexFunction::setImag(mrcpp::FunctionTree<3> *tree) {
    if (isShared()) MSG_FATAL("Cannot set in shared function");
    this->im = tree;
}

void ComplexFunction::alloc(int type) {
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (hasReal()) MSG_FATAL("Real part already allocated");
        this->re = new mrcpp::FunctionTree<3>(*MRA, this->shared_mem);
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (hasImag()) MSG_FATAL("Imaginary part already allocated");
        this->im = new mrcpp::FunctionTree<3>(*MRA, this->shared_mem);
    }
}

void ComplexFunction::free(int type) {
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (hasReal()) delete this->re;
        this->re = nullptr;
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (hasImag()) delete this->im;
        this->im = nullptr;
    }
}

/** @brief Returns the orbital meta data
     *
     * Tree sizes (nChunks) are flushed before return.
     */
FunctionData &ComplexFunction::getFunctionData() {
    this->func_data.real_size = 0;
    this->func_data.imag_size = 0;
    if (this->hasReal()) this->func_data.real_size = this->real().getNChunksUsed();
    if (this->hasImag()) this->func_data.imag_size = this->imag().getNChunksUsed();
    return this->func_data;
}

int ComplexFunction::getNNodes(int type) const {
    int nNodes = 0;
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (this->hasReal()) nNodes += this->real().getNNodes();
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (this->hasImag()) nNodes += this->imag().getNNodes();
    }
    return nNodes;
}

} // namespace mrchem
