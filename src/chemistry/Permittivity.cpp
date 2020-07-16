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

#include "chemistry/Permittivity.h"
#include "chemistry/Cavity.h"
#include <MRCPP/MWFunctions>

namespace mrchem {

Permittivity::Permittivity(const mrchem::Cavity C, double epsilon_in, double epsilon_out)
        : eps_in(epsilon_in)
        , eps_out(epsilon_out)
        , Cav(C) {}

double Permittivity::evalf(const mrcpp::Coord<3> &r) const {
    auto epsilon = eps_in * std::exp(std::log(eps_out / eps_in) * (1 - this->Cav.evalf(r)));
    if (flipped) {
        return 1 / epsilon;
    } else {
        return epsilon;
    }
}

// Not supported in C:
// C
// C_i
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// Derivative
// (1.0/4.0)*K*(K*pow(1 - C(O.x, O.y, O.z), 2)*pow(Sum(Derivative(C_i(O.x, O.y, O.z, i), O.x)/(1 - C_i(O.x, O.y, O.z,
// i)), (i, 1, N)), 2) + K*pow(1 - C(O.x, O.y, O.z), 2)*pow(Sum(Derivative(C_i(O.x, O.y, O.z, i), O.y)/(1 - C_i(O.x,
// O.y, O.z, i)), (i, 1, N)), 2) + K*pow(1 - C(O.x, O.y, O.z), 2)*pow(Sum(Derivative(C_i(O.x, O.y, O.z, i), O.z)/(1 -
// C_i(O.x, O.y, O.z, i)), (i, 1, N)), 2) + 2*(1 - C(O.x, O.y, O.z))*pow(Sum(Derivative(C_i(O.x, O.y, O.z, i), O.x)/(1 -
// C_i(O.x, O.y, O.z, i)), (i, 1, N)), 2) + 2*(1 - C(O.x, O.y, O.z))*pow(Sum(Derivative(C_i(O.x, O.y, O.z, i), O.y)/(1 -
// C_i(O.x, O.y, O.z, i)), (i, 1, N)), 2) + 2*(1 - C(O.x, O.y, O.z))*pow(Sum(Derivative(C_i(O.x, O.y, O.z, i), O.z)/(1 -
// C_i(O.x, O.y, O.z, i)), (i, 1, N)), 2) - 2*(1 - C(O.x, O.y, O.z))*Sum(Derivative(C_i(O.x, O.y, O.z, i), (O.x, 2))/(1
// - C_i(O.x, O.y, O.z, i)) + pow(Derivative(C_i(O.x, O.y, O.z, i), O.x), 2)/pow(1 - C_i(O.x, O.y, O.z, i), 2), (i, 1,
// N)) - 2*(1 - C(O.x, O.y, O.z))*Sum(Derivative(C_i(O.x, O.y, O.z, i), (O.y, 2))/(1 - C_i(O.x, O.y, O.z, i)) +
// pow(Derivative(C_i(O.x, O.y, O.z, i), O.y), 2)/pow(1 - C_i(O.x, O.y, O.z, i), 2), (i, 1, N)) - 2*(1 - C(O.x, O.y,
// O.z))*Sum(Derivative(C_i(O.x, O.y, O.z, i), (O.z, 2))/(1 - C_i(O.x, O.y, O.z, i)) + pow(Derivative(C_i(O.x, O.y, O.z,
// i), O.z), 2)/pow(1 - C_i(O.x, O.y, O.z, i), 2), (i, 1, N)))
} // namespace mrchem
