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

#include "tensor/RankOneOperator.h"

#include "MomentumOperator.h"
#include "PositionOperator.h"

namespace mrchem {

class AngularMomentumOperator final : public RankOneOperator<3> {
public:
    AngularMomentumOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D, const mrcpp::Coord<3> &o)
            : AngularMomentumOperator(PositionOperator(o), MomentumOperator(D)) {}

    AngularMomentumOperator(PositionOperator r, MomentumOperator p) {
        // Invoke operator= to assign *this operator
        RankOneOperator<3> &h = (*this);
        h[0] = (r[1] * p[2] - r[2] * p[1]);
        h[1] = (r[2] * p[0] - r[0] * p[2]);
        h[2] = (r[0] * p[1] - r[1] * p[0]);
        h[0].name() = "l[x]";
        h[1].name() = "l[y]";
        h[2].name() = "l[z]";
    }
};

} // namespace mrchem
