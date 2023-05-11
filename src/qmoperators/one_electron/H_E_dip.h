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

#include "PositionOperator.h"

/** Class H_E_dip
 *
 * @brief Electric dipole operator
 *
 * mu_i = -r_i
 *
 * The operator is based on the PositionOperator. The sign is coherent
 * with nuclear and electronic charges.
 *
 */

namespace mrchem {

class H_E_dip final : public RankOneOperator<3> {
public:
    explicit H_E_dip(const mrcpp::Coord<3> &o)
            : H_E_dip(PositionOperator(o)) {}

    explicit H_E_dip(PositionOperator r) {
        // Invoke operator= to assign *this operator
        RankOneOperator<3> &h = (*this);
        h[0] = -1.0 * r[0];
        h[1] = -1.0 * r[1];
        h[2] = -1.0 * r[2];
        h[0].name() = "h_E_dip[x]";
        h[1].name() = "h_E_dip[y]";
        h[2].name() = "h_E_dip[z]";
    }
};

} // namespace mrchem
