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
#include "qmfunction_fwd.h"

namespace mrchem {

/* POD struct for function meta data. Used for simple MPI communication. */
struct FunctionData {
    int nChunksReal;
    int nChunksImag;
    bool conjugate;
};

class QMFunction {
public:
    QMFunction(bool share, mrcpp::FunctionTree<3> *r = nullptr, mrcpp::FunctionTree<3> *i = nullptr);
    QMFunction(const QMFunction &func);
    QMFunction &operator=(const QMFunction &func);
    virtual ~QMFunction();

    void set(int type, mrcpp::FunctionTree<3> *func);
    void alloc(int type);
    void free(int type = NUMBER::Total);

    int getNNodes(int type = NUMBER::Total) const;
    bool conjugate() const { return this->func_data.conjugate; }
    FunctionData &getFunctionData();

    double norm() const;
    double squaredNorm() const;

    void add(ComplexDouble c, QMFunction inp);
    void rescale(ComplexDouble c);
    void crop(double prec);

    bool isShared() const { return (this->shared_mem == nullptr) ? false : true; }
    bool hasReal() const { return (this->re == nullptr) ? false : true; }
    bool hasImag() const { return (this->im == nullptr) ? false : true; }

    mrcpp::FunctionTree<3> &real() { return *this->re; }
    mrcpp::FunctionTree<3> &imag() { return *this->im; }

    const mrcpp::FunctionTree<3> &real() const { return *this->re; }
    const mrcpp::FunctionTree<3> &imag() const { return *this->im; }

protected:
    FunctionData func_data;
    mrcpp::SharedMemory *shared_mem;
    mrcpp::FunctionTree<3> *re; ///< Real part of function
    mrcpp::FunctionTree<3> *im; ///< Imaginary part of function
};

} // namespace mrchem
