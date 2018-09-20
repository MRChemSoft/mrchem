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

class GeometryDerivatives final {
public:
    GeometryDerivatives(int k) {
        this->nuclear = DoubleMatrix::Zero(k, 3);
        this->electronic = DoubleMatrix::Zero(k, 3);
    }
    ~GeometryDerivatives() { }

    DoubleMatrix get() const { return this->nuclear + this->electronic; }
    DoubleMatrix &getNuclear() { return this->nuclear; }
    DoubleMatrix &getElectronic() { return this->electronic; }

    friend std::ostream& operator<<(std::ostream &o, const GeometryDerivatives &geomderiv) {
        DoubleMatrix gd = geomderiv.get();

        int oldPrec = mrcpp::Printer::setPrecision(10);
        double au = gd.norm();
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                   Geometry derivatives                     "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Length of vector:   (au)    " << std::setw(30) << au        <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        for (int k = 0; k < gd.rows(); k++) {
            o<<std::setw(19)<<gd(k,0)<<std::setw(20)<<gd(k,1)<<std::setw(20)<<gd(k,2)<<std::endl;
        }
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                                                            "<<std::endl;
        mrcpp::Printer::setPrecision(oldPrec);
        return o;
    }
protected:
    DoubleMatrix nuclear;
    DoubleMatrix electronic;
};

} //namespace mrchem
