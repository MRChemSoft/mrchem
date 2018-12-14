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
} // namespace mrchem
