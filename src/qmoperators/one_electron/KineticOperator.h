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

#include "tensor/RankZeroOperator.h"

#include "MomentumOperator.h"

/** @class KineticOperator
 *
 * @brief Operator for kinetic energy
 *
 * This operator is constructed as the square of the more fundamental
 * MomentumOperator. The general base class functions for calculation of
 * expectation values are overwritten, as they can be improved due to
 * symmetry.
 *
 */

namespace mrchem {

class KineticOperator final : public RankZeroOperator {
public:
    explicit KineticOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D)
            : KineticOperator(MomentumOperator(D)) {}

    explicit KineticOperator(MomentumOperator p) {
        // Invoke operator= to assign *this operator
        RankZeroOperator &t = (*this);
        t = 0.5 * (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        t.name() = "T";
    }
};

} // namespace mrchem
