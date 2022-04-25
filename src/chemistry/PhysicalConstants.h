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
#include "MRCPP/Printer"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace mrchem {

class PhysicalConstants {
public:
    static PhysicalConstants &Initialize(const json &constants);
    static double get(const std::string &key) {
        try {
            if (initialized) {
                return constants_[key];
            } else {
                return testConstants[key];
            }
        } catch (...) { MSG_ABORT("Error getting constant with name: " + key); }
    }

    static void Print();

    PhysicalConstants() = default;
    ~PhysicalConstants() = default;

    PhysicalConstants(const PhysicalConstants &) = delete;
    PhysicalConstants &operator=(const PhysicalConstants &) = delete;
    PhysicalConstants &operator=(const PhysicalConstants &&) = delete;
    PhysicalConstants(PhysicalConstants &&) = delete;

    static bool initialized;

private:
    PhysicalConstants(const json &constants) { constants_ = constants; }
    static json constants_;
    static json testConstants;
};

} // namespace mrchem
