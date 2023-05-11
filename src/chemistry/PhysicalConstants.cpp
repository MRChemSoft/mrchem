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

#include "PhysicalConstants.h"
#include "MRCPP/Printer"
#include "utils/print_utils.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace mrchem {

bool PhysicalConstants::initialized = false;
json PhysicalConstants::constants_ = json();

// clang-format off
json PhysicalConstants::testConstants = {
    {"angstrom2bohrs", 1.8897261246257702},
    {"dipmom_au2debye", 2.5417464739297717},
    {"electron_g_factor", -2.00231930436256},
    {"fine_structure_constant", 0.0072973525693},
    {"hartree2ev", 27.211386245988},
    {"hartree2kcalmol", 627.5094740630558},
    {"hartree2kjmol", 2625.4996394798254},
    {"hartree2simagnetizability", 78.9451185},
    {"hartree2wavenumbers", 219474.6313632},
    {"light_speed", 137.035999084}
};
// clang-format on

PhysicalConstants &PhysicalConstants::Initialize(const json &constants) {
    initialized = true;
    static PhysicalConstants obj(constants);
    return obj;
}

/** @brief Pretty print physical constants */
void PhysicalConstants::Print() {
    mrcpp::print::separator(0, '~');
    print_utils::centered_text(0, "Physical Constants");
    mrcpp::print::separator(0, '~');
    print_utils::json(0, constants_, true);
    mrcpp::print::separator(0, '~', 2);
}

} // namespace mrchem
