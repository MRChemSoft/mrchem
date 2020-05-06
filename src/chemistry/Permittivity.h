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

#include "chemistry/Cavity.h"
#include <MRCPP/MWFunctions>

namespace mrchem {
class Cavity;

class Permittivity final : public mrcpp::RepresentableFunction<3> {
public:
    Permittivity(const Cavity C, double epsilon_in, double epsilon_out = 1.0);
    double evalf(const mrcpp::Coord<3> &r) const override;

    void flipFunction(bool is_flipped) { this->flipped = is_flipped; };
    bool isInverse() { return this->flipped; }
    std::vector<mrcpp::Coord<3>> getCoordinates() const { return this->Cav.getCoordinates(); }
    std::vector<double> getRadius() const { return this->Cav.getRadius(); }
    auto getGradVector() const { return this->Cav.getGradVector(); }
    double eps_in;
    double eps_out;

private:
    bool flipped = false;
    Cavity Cav;
};

} // namespace mrchem
