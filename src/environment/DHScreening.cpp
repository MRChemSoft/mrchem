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

#include "DHScreening.h"

namespace mrchem {

DHScreening::DHScreening(std::shared_ptr<mrchem::Cavity> cavity_ion, double kappa_out, const std::string &formulation)
        : StepFunction(cavity_ion, 0.0, kappa_out)
        , formulation(formulation) {}

double DHScreening::evalf(const mrcpp::Coord<3> &r) const {
    auto c_pin = this->cavity_p;
    auto kappa_squared = std::pow(this->out, 2) * (1 - c_pin->evalf(r));
    return kappa_squared;
}

} // namespace mrchem
