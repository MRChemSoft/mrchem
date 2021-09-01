/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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
#include "ZoraOperator.h"

namespace mrchem {

class ZoraKineticOperator final : public RankZeroOperator {
public:
    ZoraKineticOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D, ZoraOperator kappa)
            : ZoraKineticOperator(MomentumOperator(D), kappa) {}

    ZoraKineticOperator(MomentumOperator p, ZoraOperator kappa) {
        // Invoke operator= to assign *this operator
        RankZeroOperator &t = (*this);
        t = 0.5 * (p[0] * kappa * p[0] + p[1] * kappa * p[1] + p[2] * kappa * p[2]);
        t.name() = "T_zora";
    }
};

} // namespace mrchem
