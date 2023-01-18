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

#include "qmoperators/QMDerivative.h"

namespace mrchem {

class MomentumOperator final : public RankOneOperator<3> {
public:
    MomentumOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D) {
        auto p_x = std::make_shared<QMDerivative>(0, D, true);
        auto p_y = std::make_shared<QMDerivative>(1, D, true);
        auto p_z = std::make_shared<QMDerivative>(2, D, true);

        // Invoke operator= to assign *this operator
        RankOneOperator<3> &p = (*this);
        p[0] = p_x;
        p[1] = p_y;
        p[2] = p_z;
        p[0].name() = "p[x]";
        p[1].name() = "p[y]";
        p[2].name() = "p[z]";
    }
};

} // namespace mrchem
