/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "NuclearFunction.h"

#include "utils/math_utils.h"

namespace mrchem {

class FiniteNucleusSphere : public NuclearFunction {
public:
    FiniteNucleusSphere() = default;

    double evalf(const mrcpp::Coord<3> &r) const override {
        double result = 0.0;
        for (int i = 0; i < this->nuclei.size(); i++) {
            const auto &R = this->nuclei[i].getCoord();
            auto R1 = math_utils::calc_distance(r, R);
            auto Z = this->nuclei[i].getCharge();
            auto R0 = this->param2[i];
            if (R1 <= R0) {
                result += -(Z / (2.0*R0)) * (3.0 - (R1*R1)/(R0*R0));
            } else { 
                result += -(Z / R1);
            }
        }
        return result;
    }

    std::string getParamName1() const override{ return "RMS"; }
    std::string getParamName2() const override{ return "R0"; }
    double calcParam1(double prec, const Nucleus &nuc) const override { return nuc.getRMSRadius(); }
    double calcParam2(double prec, const Nucleus &nuc) const override {
        auto RMS = nuc.getRMSRadius();
        auto RMS2 = RMS*RMS;
        auto R0 = std::sqrt(RMS2*(5.0/3.0));
        return R0;
    }
};

} // namespace mrchem
