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

#pragma once

#include "mrchem.h"
#include "parallel.h"

namespace mrchem {

struct FunctionData {
    int real_size;
    int imag_size;
    bool is_shared;
};

class ComplexFunction final {
public:
    explicit ComplexFunction(bool share);
    ~ComplexFunction();

    void alloc(int type);
    void free(int type);

    bool isShared() const { return this->func_data.is_shared; }
    bool hasReal() const { return (this->re == nullptr) ? false : true; }
    bool hasImag() const { return (this->im == nullptr) ? false : true; }

    int getNNodes(int type) const;
    FunctionData &getFunctionData();

    void setReal(mrcpp::FunctionTree<3> *tree);
    void setImag(mrcpp::FunctionTree<3> *tree);

    mrcpp::FunctionTree<3> &real() { return *this->re; }
    mrcpp::FunctionTree<3> &imag() { return *this->im; }

    const mrcpp::FunctionTree<3> &real() const { return *this->re; }
    const mrcpp::FunctionTree<3> &imag() const { return *this->im; }

private:
    FunctionData func_data;
    mrcpp::SharedMemory *shared_mem;
    mrcpp::FunctionTree<3> *re; ///< Real part of function
    mrcpp::FunctionTree<3> *im; ///< Imaginary part of function
};

} // namespace mrchem
