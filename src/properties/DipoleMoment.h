/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "chemistry/PhysicalConstants.h"

#include "utils/math_utils.h"
#include "utils/print_utils.h"

namespace mrchem {

// clang-format off
class DipoleMoment final {
public:
    explicit DipoleMoment(const mrcpp::Coord<3> &o = {}) : r_O(o) {}

    void setOrigin(const mrcpp::Coord<3> &o) { this->r_O = o; }
    const mrcpp::Coord<3> &getOrigin() const { return this->r_O; }

    DoubleVector getTensor() const { return getNuclear() + getElectronic(); }
    DoubleVector &getNuclear() { return this->nuc_tensor; }
    DoubleVector &getElectronic() { return this->el_tensor; }
    const DoubleVector &getNuclear() const { return this->nuc_tensor; }
    const DoubleVector &getElectronic() const { return this->el_tensor; }

    void print(const std::string &id) const {
        auto el_au = getElectronic().norm();
        auto nuc_au = getNuclear().norm();
        auto tot_au = getTensor().norm();

        auto el_db = el_au * PhysicalConstants::get("dipmom_au2debye");
        auto nuc_db = nuc_au * PhysicalConstants::get("dipmom_au2debye");
        auto tot_db = tot_au * PhysicalConstants::get("dipmom_au2debye");

        mrcpp::print::header(0, "Dipole Moment (" + id + ")");
        print_utils::coord(0, "r_O", getOrigin());
        mrcpp::print::separator(0, '-');
        print_utils::vector(0, "Electronic vector", getElectronic());
        print_utils::scalar(0, "Magnitude", el_au, "(au)");
        print_utils::scalar(0, "         ", el_db, "(Debye)");
        mrcpp::print::separator(0, '-');
        print_utils::vector(0, "Nuclear vector", getNuclear());
        print_utils::scalar(0, "Magnitude", nuc_au, "(au)");
        print_utils::scalar(0, "         ", nuc_db, "(Debye)");
        mrcpp::print::separator(0, '-');
        print_utils::vector(0, "Total vector", getTensor());
        print_utils::scalar(0, "Magnitude", tot_au, "(au)");
        print_utils::scalar(0, "         ", tot_db, "(Debye)");
        mrcpp::print::separator(0, '=', 2);
    }

    nlohmann::json json() const {
        return {
            {"r_O", getOrigin()},
            {"vector_nuc", print_utils::eigen_to_vector(getNuclear(), 1.0e-12)},
            {"vector_el", print_utils::eigen_to_vector(getElectronic(), 1.0e-12)},
            {"vector", print_utils::eigen_to_vector(getTensor(), 1.0e-12)},
            {"magnitude", getTensor().norm() }
        };
    }

private:
    mrcpp::Coord<3> r_O;
    DoubleVector nuc_tensor{math_utils::init_nan(3)};
    DoubleVector el_tensor{math_utils::init_nan(3)};
};
// clang-format on

} // namespace mrchem
