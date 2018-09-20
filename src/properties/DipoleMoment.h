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

namespace mrchem {

class DipoleMoment final {
public:
    DipoleMoment() {
        this->nuclear = DoubleVector::Zero(3);
        this->electronic = DoubleVector::Zero(3);
    }
    ~DipoleMoment() { }

    DoubleVector get() const { return this->nuclear + this->electronic; }
    DoubleVector &getNuclear() { return this->nuclear; }
    DoubleVector &getElectronic() { return this->electronic; }

    friend std::ostream& operator<<(std::ostream &o, const DipoleMoment &dipole) {
        DoubleVector mu = dipole.get();

        int oldPrec = mrcpp::Printer::setPrecision(10);
        double au = mu.norm();
        double debye = au * PHYSCONST::Debye;
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                         Dipole moment                      "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Length of vector:   (au)    " << std::setw(30) << au        <<std::endl;
        o<<"                     (Debye) " << std::setw(30) << debye     <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<std::setw(19)<<mu(0)<<std::setw(20)<<mu(1)<<std::setw(20)<<mu(2)<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                                                            "<<std::endl;
        mrcpp::Printer::setPrecision(oldPrec);
        return o;
    }
protected:
    DoubleVector nuclear;
    DoubleVector electronic;
};

} //namespace mrchem
