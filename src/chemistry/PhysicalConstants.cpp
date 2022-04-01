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

#include "PhysicalConstants.h"

#include <nlohmann/json.hpp>

namespace mrchem {

PhysicalConstants *PhysicalConstants::_instance = NULL;

void PhysicalConstants::setConstants(const nlohmann::json &constants) {
    hasData = true;
    _constants = constants;
}
double PhysicalConstants::get(std::string name) {
    if (hasData) {
        return _constants[name];
    } else {
        return testConstants[name];
    }
}

PhysicalConstants::PhysicalConstants() {
    testConstants["pi"] = 3.141592653589793;
    testConstants["pi_sqrt"] = 1.772453850905516;
    testConstants["electron_g_factor"] = -2.00231930436256;
    testConstants["fine_structure_constant"] = 0.0072973525693;
    testConstants["hartree2kjmol"] = 2625.4996394798254;
    testConstants["hartree2kcalmol"] = 627.5094740630558;
    testConstants["hartree2ev"] = 27.211386245988;
    testConstants["hartree2simagnetizability"] = 78.9451185;
}

PhysicalConstants *PhysicalConstants::getInstance() {
    if (_instance == NULL) { _instance = new PhysicalConstants(); }
    return _instance;
}

} // namespace mrchem
