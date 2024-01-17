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

#include "qmoperators/one_electron/NablaOperator.h"

namespace mrchem {

class MomentumOperator final : public RankOneOperator<3> {
public:
    MomentumOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D)
            : MomentumOperator(NablaOperator(D, true)) {}

    MomentumOperator(NablaOperator D) {
        // Invoke operator= to assign *this operator
        RankOneOperator<3> &p = (*this);
        p[0] = -1.0*D[0];
        p[1] = -1.0*D[1];
        p[2] = -1.0*D[2];
        p[0].name() = "p[x]";
        p[1].name() = "p[y]";
        p[2].name() = "p[z]";
    }
};

} // namespace mrchem
