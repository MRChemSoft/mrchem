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
 
#include "mrchem.h"
#include "utils/math_utils.h"
#include "utils/print_utils.h"

#include <nlohmann/json.hpp>
#include <string>

namespace mrchem {

class HirshfeldCharges final {
public:
    DoubleVector getVector() const { return this->hirshfeld_charges; }

    void setVector(const DoubleVector &v) { this->hirshfeld_charges = v; }

    void print(const std::string &id) const {
        mrcpp::print::header(0, "Hirshfeld Charges (" + id + ")");
        mrcpp::print::separator(0, '-');
        for (int i = 0; i < hirshfeld_charges.size(); i++) {
            std::string text = "Charge of atom " + std::to_string(i);
            print_utils::scalar(0, text, hirshfeld_charges(i));
        }
        mrcpp::print::separator(0, '-');
        print_utils::scalar(0, "Sum of Hirshfeld charges", getVector().sum(), "(au)", -1, true);
        mrcpp::print::separator(0, '=');
    }

    nlohmann::json json() const { return {{"total", print_utils::eigen_to_vector(getVector(), 1.0e-12)}}; }

protected:
    DoubleVector hirshfeld_charges;
};

} // namespace mrchem