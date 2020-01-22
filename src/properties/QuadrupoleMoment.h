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

#include "mrchem.h"

#include "utils/math_utils.h"
#include "utils/print_utils.h"

namespace mrchem {

// clang-format off
class QuadrupoleMoment final {
public:
    DoubleMatrix getTensor() const { return getNuclear() + getElectronic(); }
    DoubleMatrix &getNuclear() { return this->nuc_tensor; }
    DoubleMatrix &getElectronic() { return this->el_tensor; }
    const DoubleMatrix &getNuclear() const { return this->nuc_tensor; }
    const DoubleMatrix &getElectronic() const { return this->el_tensor; }

    void print() const {
        mrcpp::print::header(0, "Quadrupole Moment");
        print_utils::matrix(0, "Electronic tensor", getElectronic());
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Nuclear tensor", getNuclear());
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Total tensor", getTensor());
        mrcpp::print::separator(0, '=', 2);
    }

    nlohmann::json json() const {
        return {
            {"tensor_nuc", math_utils::eigen_to_vector(getNuclear())},
            {"tensor_el", math_utils::eigen_to_vector(getElectronic())},
            {"tensor", math_utils::eigen_to_vector(getTensor())}
        };
    }

private:
    DoubleMatrix nuc_tensor{math_utils::init_nan(3)};
    DoubleMatrix el_tensor{math_utils::init_nan(3)};
};
// clang-format on

} // namespace mrchem
