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

class PointNucleusHFYGB : public NuclearFunction {
public:
    PointNucleusHFYGB() = default;

    double evalf(const mrcpp::Coord<3> &r) const override {
        double result = 0.0;
        for (int i = 0; i < this->nuclei.size(); i++) {
            const auto &R = this->nuclei[i].getCoord();
            double R1 = math_utils::calc_distance(R, r);
            double Z = this->nuclei[i].getCharge();
            double S_i = this->param2[i];
            R1 /= S_i;
            double c = -1.0 / (3.0 * mrcpp::root_pi);
            double partResult = -std::erf(R1) / R1 + c * (std::exp(-R1 * R1) + 16.0 * std::exp(-4.0 * R1 * R1));
            result += Z * partResult / S_i;
        }
        return result;
    }

    std::string getParamName1() const override { return "Precision"; }
    std::string getParamName2() const override { return "Smoothing"; }
    double calcParam1(double prec, const Nucleus &nuc) const override { return prec; }
    double calcParam2(double prec, const Nucleus &nuc) const override {
        auto Z = nuc.getCharge();
        double tmp = 0.00435 * prec / std::pow(Z, 5.0);
        return std::cbrt(tmp);
    }
};

} // namespace mrchem
