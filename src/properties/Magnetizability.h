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

class Magnetizability final {
public:
    Magnetizability() {
        this->diamagnetic = DoubleMatrix::Zero(3,3);
        this->paramagnetic = DoubleMatrix::Zero(3,3);
    }
    ~Magnetizability() { }

    DoubleMatrix get() const { return this->diamagnetic + this->paramagnetic; }
    DoubleMatrix &getDiamagnetic() { return this->diamagnetic; }
    DoubleMatrix &getParamagnetic() { return this->paramagnetic; }

    friend std::ostream& operator<<(std::ostream &o, const Magnetizability &mag) {
        DoubleMatrix dia = mag.diamagnetic;
        DoubleMatrix para = mag.paramagnetic;
        DoubleMatrix tot = dia + para;

        double w_au = 0.0;  // Only static magnetizability
        double isoDMau = dia.trace()/3.0;
        double isoPMau = para.trace()/3.0;
        double isoTMau = isoDMau + isoPMau;

        double isoDMsi = isoDMau * PHYSCONST::JT_m2; // SI units (J/T^2 10^{-30})
        double isoPMsi = isoPMau * PHYSCONST::JT_m2; // SI units (J/T^2 10^{-30})
        double isoTMsi = isoTMau * PHYSCONST::JT_m2; // SI units (J/T^2 10^{-30})

        int oldPrec = mrcpp::Printer::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                   Magnetizability tensor                   "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Frequency:       (au)       " << std::setw(30) << w_au      <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"-------------------- Isotropic averages --------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Total            (au)       " << std::setw(30) << isoTMau   <<std::endl;
        o<<" Diamagnetic      (au)       " << std::setw(30) << isoDMau   <<std::endl;
        o<<" Paramagnetic     (au)       " << std::setw(30) << isoPMau   <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Total            (SI)       " << std::setw(30) << isoTMsi   <<std::endl;
        o<<" Diamagnetic      (SI)       " << std::setw(30) << isoDMsi   <<std::endl;
        o<<" Paramagnetic     (SI)       " << std::setw(30) << isoPMsi   <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"-------------------------- Total ---------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<std::setw(19)<<tot(0,0)<<std::setw(20)<<tot(1,0)<<std::setw(20)<<tot(2,0)<<std::endl;
        o<<std::setw(19)<<tot(1,0)<<std::setw(20)<<tot(1,1)<<std::setw(20)<<tot(2,1)<<std::endl;
        o<<std::setw(19)<<tot(2,0)<<std::setw(20)<<tot(1,2)<<std::setw(20)<<tot(2,2)<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"----------------------- Diamagnetic ------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<std::setw(19)<<dia(0,0)<<std::setw(20)<<dia(1,0)<<std::setw(20)<<dia(2,0)<<std::endl;
        o<<std::setw(19)<<dia(1,0)<<std::setw(20)<<dia(1,1)<<std::setw(20)<<dia(2,1)<<std::endl;
        o<<std::setw(19)<<dia(2,0)<<std::setw(20)<<dia(1,2)<<std::setw(20)<<dia(2,2)<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"----------------------- Paramagnetic -----------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<std::setw(19)<<para(0,0)<<std::setw(20)<<para(1,0)<<std::setw(20)<<para(2,0)<<std::endl;
        o<<std::setw(19)<<para(1,0)<<std::setw(20)<<para(1,1)<<std::setw(20)<<para(2,1)<<std::endl;
        o<<std::setw(19)<<para(2,0)<<std::setw(20)<<para(1,2)<<std::setw(20)<<para(2,2)<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                                                            "<<std::endl;
        mrcpp::Printer::setPrecision(oldPrec);
        return o;
    }
protected:
    DoubleMatrix diamagnetic;
    DoubleMatrix paramagnetic;
    static const std::string name;
};

} //namespace mrchem
