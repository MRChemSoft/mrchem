/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

/** @class SCFEnergy
 *
 * @brief Simple POD container to hold the different contributions to the SCF energy.
 *
 */

namespace mrchem {

// clang-format off
class SCFEnergy final {
public:
    SCFEnergy(double nuc = 0.0, double el = 0.0, double orb = 0.0,
              double kin = 0.0, double en = 0.0, double ee = 0.0,
              double xc = 0.0, double x = 0.0, double nex = 0.0,
              double ext = 0.0, double r = 0.0) :
        E_nuc(nuc), E_el(el), E_orb(orb), E_kin(kin), E_en(en),
        E_ee(ee), E_x(x), E_xc(xc), E_nex(nex), E_ext(ext), E_r(r){ }

    double getTotalEnergy() const { return this->E_nuc + this->E_el; }
    double getNuclearEnergy() const { return this->E_nuc; }
    double getElectronicEnergy() const { return this->E_el; }

    double getOrbitalEnergy() const { return this->E_orb; }
    double getKineticEnergy() const { return this->E_kin; }
    double getElectronNuclearEnergy() const { return this->E_en; }
    double getElectronElectronEnergy() const { return this->E_ee; }
    double getExchangeCorrelationEnergy() const { return this->E_xc; }
    double getExchangeEnergy() const { return this->E_x; }
    double getReactionEnergy() const {return this->E_r; }

    friend std::ostream& operator<<(std::ostream &o, const SCFEnergy &en) {
        auto E_au   = en.E_nuc + en.E_el;
        auto E_kJ   = E_au * PHYSCONST::kJ;
        auto E_kcal = E_au * PHYSCONST::kcal;
        auto E_eV   = E_au * PHYSCONST::eV;

        auto oldPrec = mrcpp::Printer::setPrecision(15);
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                         SCF Energy                         " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o << " Sum orbital energy:          " << std::setw(29) << en.E_orb  << std::endl;
        o << " Kinetic energy:              " << std::setw(29) << en.E_kin  << std::endl;
        o << " E-N energy:                  " << std::setw(29) << en.E_en   << std::endl;
        o << " Coulomb energy:              " << std::setw(29) << en.E_ee   << std::endl;
        o << " Exchange energy:             " << std::setw(29) << en.E_x    << std::endl;
        o << " X-C energy:                  " << std::setw(29) << en.E_xc   << std::endl;
        o << " El. External field energy:   " << std::setw(29) << en.E_ext  << std::endl;
        o << " Nuc. External field energy:  " << std::setw(29) << en.E_nex  << std::endl;
        o << " Reac. field energy:          " << std::setw(29) << en.E_r    << std::endl;
        o << "                                                            " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o << " Electronic energy:           " << std::setw(29) << en.E_el   << std::endl;
        o << " Nuclear energy:              " << std::setw(29) << en.E_nuc  << std::endl;
        o << "                                                            " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o << " Total energy       (au)      " << std::setw(29) << E_au      << std::endl;
        o << "                    (kJ/mol)  " << std::setw(29) << E_kJ      << std::endl;
        o << "                    (kcal/mol)" << std::setw(29) << E_kcal    << std::endl;
        o << "                    (eV)      " << std::setw(29) << E_eV      << std::endl;
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                                                            " << std::endl;
        mrcpp::Printer::setPrecision(oldPrec);

        return o;
    }

private:
    double E_nuc{0.0};
    double E_el{0.0};

    double E_orb;
    double E_kin;
    double E_en;
    double E_ee;
    double E_x;
    double E_xc;
    double E_nex;
    double E_ext;
    double E_r;
};
// clang-format on

} // namespace mrchem
