/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include <nlohmann/json.hpp>

#include "mrchem.h"
#include "utils/print_utils.h"

/** @class SCFEnergy
 *
 * @brief Simple POD container to hold the different contributions to the SCF energy.
 *
 */

namespace mrchem {

// clang-format off
class SCFEnergy final {
public:
    explicit SCFEnergy(double kin = 0.0, double nn = 0.0,
                       double en = 0.0, double ee = 0.0,
                       double x = 0.0, double xc = 0.0,
                       double next = 0.0, double eext = 0.0) :
        E_kin(kin), E_nn(nn), E_en(en), E_ee(ee),
        E_x(x), E_xc(xc), E_next(next), E_eext(eext) {
            E_nuc = E_nn + E_next;
            E_el = E_kin + E_en + E_ee + E_xc + E_x + E_eext;
        }

    double getTotalEnergy() const { return this->E_nuc + this->E_el; }
    double getNuclearEnergy() const { return this->E_nuc; }
    double getElectronicEnergy() const { return this->E_el; }

    double getKineticEnergy() const { return this->E_kin; }
    double getNuclearNuclearEnergy() const { return this->E_nn; }
    double getElectronNuclearEnergy() const { return this->E_en; }
    double getElectronElectronEnergy() const { return this->E_ee; }
    double getElectronExternalEnergy() const { return this->E_eext; }
    double getNuclearExternalEnergy() const { return this->E_next; }
    double getExchangeCorrelationEnergy() const { return this->E_xc; }
    double getExchangeEnergy() const { return this->E_x; }

    void print(const std::string &id) const {
        auto E_au = E_nuc + E_el;
        auto E_eV = E_au * PHYSCONST::eV;
        auto E_kJ = E_au * PHYSCONST::kJ;
        auto E_kcal = E_au * PHYSCONST::kcal;

        auto pprec = 2 * mrcpp::Printer::getPrecision();
        mrcpp::print::header(0, "Molecular Energy (" + id + ")");
        print_utils::scalar(0, "Kinetic energy   ", E_kin,  "(au)", pprec, false);
        print_utils::scalar(0, "E-N energy       ", E_en,   "(au)", pprec, false);
        print_utils::scalar(0, "Coulomb energy   ", E_ee,   "(au)", pprec, false);
        print_utils::scalar(0, "Exchange energy  ", E_x,    "(au)", pprec, false);
        print_utils::scalar(0, "X-C energy       ", E_xc,   "(au)", pprec, false);
        print_utils::scalar(0, "Ext. field (el)  ", E_eext, "(au)", pprec, false);
        mrcpp::print::separator(0, '-');
        print_utils::scalar(0, "N-N energy       ", E_nn,   "(au)", pprec, false);
        print_utils::scalar(0, "Ext. field (nuc) ", E_next, "(au)", pprec, false);
        mrcpp::print::separator(0, '-');
        print_utils::scalar(0, "Electronic energy", E_el,   "(au)", pprec, false);
        print_utils::scalar(0, "Nuclear energy   ", E_nuc,  "(au)", pprec, false);
        mrcpp::print::separator(0, '-');
        print_utils::scalar(0, "Total energy     ", E_au,   "(au)", pprec, true);
        print_utils::scalar(0, "                 ", E_kcal, "(kcal/mol)", pprec, true);
        print_utils::scalar(0, "                 ", E_kJ,   "(kJ/mol)", pprec, true);
        print_utils::scalar(0, "                 ", E_eV,   "(eV)", pprec, true);
        mrcpp::print::separator(0, '=', 2);
    }

    nlohmann::json json() const {
        return {
            {"E_kin", E_kin},
            {"E_nn", E_nn},
            {"E_en", E_en},
            {"E_ee", E_ee},
            {"E_next", E_next},
            {"E_eext", E_eext},
            {"E_x", E_x},
            {"E_xc", E_xc},
            {"E_el", E_el},
            {"E_nuc", E_nuc},
            {"E_tot", E_nuc + E_el}
        };
    }

private:
    double E_nuc{0.0};
    double E_el{0.0};

    double E_kin{0.0};
    double E_nn{0.0};
    double E_en{0.0};
    double E_ee{0.0};
    double E_x{0.0};
    double E_xc{0.0};
    double E_next{0.0};
    double E_eext{0.0};
};
// clang-format on

} // namespace mrchem
