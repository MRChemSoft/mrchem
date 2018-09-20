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

#include "QMFunction.h"

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

QMFunction& QMFunction::operator=(const QMFunction &func) {
    if (this != &func) {
        this->re = func.re;
        this->im = func.im;
    }
    return *this;
}

void QMFunction::alloc(int type) {
    if (type == NUMBER::Real or type == NUMBER::Total) {
        if (this->hasReal()) MSG_FATAL("Function not empty");
        this->re = new mrcpp::FunctionTree<3>(*MRA);
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        if (this->hasImag()) MSG_FATAL("Function not empty");
        this->im = new mrcpp::FunctionTree<3>(*MRA);
    }
}

void QMFunction::clear(int type) {
    if (type == NUMBER::Real or type == NUMBER::Total) {
        this->re = nullptr;
    }
    if (type == NUMBER::Imag or type == NUMBER::Total) {
        this->im = nullptr;
    }
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

} //namespace mrchem
