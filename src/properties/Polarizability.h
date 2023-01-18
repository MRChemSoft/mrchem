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

#include <nlohmann/json.hpp>

#include "chemistry/PhysicalConstants.h"
#include "mrchem.h"

#include "utils/math_utils.h"
#include "utils/print_utils.h"

namespace mrchem {

// clang-format off
class Polarizability final {
public:
    explicit Polarizability(double omega = 0.0, const mrcpp::Coord<3> &o = {}) : frequency(omega), r_O(o) {}

    void setFrequency(double omega) { this->frequency = omega; }
    void setOrigin(const mrcpp::Coord<3> &o) { this->r_O = o; }

    double getFrequency() const { return this->frequency; }
    const mrcpp::Coord<3> &getOrigin() const { return this->r_O; }

    DoubleMatrix &getTensor() { return this->tensor; }
    const DoubleMatrix &getTensor() const { return this->tensor; }

    void print(const std::string &id) const {
        auto w_au = getFrequency();
        auto w_cm = PhysicalConstants::get("hartree2wavenumbers") * w_au;
        auto dynamic = (w_au > mrcpp::MachineZero);
        auto l_nm = (dynamic) ? (1.0e7 / w_cm) : 0.0;
        auto iso_au = getTensor().trace() / 3.0;

        mrcpp::print::header(0, "Polarizability (" + id + ")");
        if (dynamic) print_utils::scalar(0, "Wavelength", l_nm, "(nm)");
        print_utils::scalar(0, "Frequency", w_au, "(au)");
        print_utils::scalar(0, "         ", w_cm, "(cm-1)");
        print_utils::coord(0, "r_O", getOrigin());
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Total tensor", getTensor());
        print_utils::scalar(0, "Isotropic average", iso_au, "(au)");
        mrcpp::print::separator(0, '=', 2);
    }

    nlohmann::json json() const {
        return {
            {"r_O", getOrigin()},
            {"frequency", getFrequency()},
            {"tensor", print_utils::eigen_to_vector(getTensor(), 1.0e-12)},
            {"isotropic_average", getTensor().trace() / 3.0 }
        };
    }

private:
    double frequency;
    mrcpp::Coord<3> r_O;
    DoubleMatrix tensor{math_utils::init_nan(3,3)};
};
// clang-format on

} // namespace mrchem
