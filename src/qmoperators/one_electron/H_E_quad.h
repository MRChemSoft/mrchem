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

#include "tensor/RankTwoOperator.h"

#include "PositionOperator.h"

/** Class H_E_quad
 *
 * @brief Electric quadrupole operator
 *
 * Q_ij = -3/2*r_i*r_j + d_ij/2*r^2
 *
 * The operator is based on the PositionOperator. The sign is coherent
 * with nuclear and electronic charges.
 *
 */
namespace mrchem {

class H_E_quad final : public RankTwoOperator<3, 3> {
public:
    explicit H_E_quad(const mrcpp::Coord<3> &o)
            : H_E_quad(PositionOperator(o)) {}

    explicit H_E_quad(PositionOperator r) {
        // Invoke operator= to assign *this operator
        RankTwoOperator &h = (*this);
        h[0][0] = -1.0 * r[0](r[0]) + 0.5 * r[1](r[1]) + 0.5 * r[2](r[2]);
        h[0][1] = -1.5 * r[0](r[1]);
        h[0][2] = -1.5 * r[0](r[2]);
        h[1][0] = -1.5 * r[1](r[0]);
        h[1][1] = -1.0 * r[1](r[1]) + 0.5 * r[0](r[0]) + 0.5 * r[2](r[2]);
        h[1][2] = -1.5 * r[1](r[2]);
        h[2][0] = -1.5 * r[2](r[0]);
        h[2][1] = -1.5 * r[2](r[1]);
        h[2][2] = -1.0 * r[2](r[2]) + 0.5 * r[0](r[0]) + 0.5 * r[1](r[1]);
        h[0][0].name() = "h_e_quad[x,x]";
        h[0][1].name() = "h_e_quad[x,y]";
        h[0][2].name() = "h_e_quad[x,z]";
        h[1][0].name() = "h_e_quad[y,x]";
        h[1][1].name() = "h_e_quad[y,y]";
        h[1][2].name() = "h_e_quad[y,z]";
        h[2][0].name() = "h_e_quad[z,x]";
        h[2][1].name() = "h_e_quad[z,y]";
        h[2][2].name() = "h_e_quad[z,z]";
    }
};

} // namespace mrchem
