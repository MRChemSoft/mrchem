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

namespace mrchem {

/** @class H_BB_dia
 *
 * @brief Diamagnetic magnetizability operator
 *
 * Interaction operator obtained by differentiating the spin Hamiltonian twice
 * wrt the external magnetic field B:
 *
 * d^2H/dB^2 = H_BB_dia
 *
 * H_BB_dia = \sum_j (r_{jO} \cdot r_{jO})1 - r_{jO}r_{jO}^T
 */

class H_BB_dia final : public RankTwoOperator<3, 3> {
public:
    explicit H_BB_dia(const mrcpp::Coord<3> &o)
            : H_BB_dia(PositionOperator(o)) {}

    explicit H_BB_dia(PositionOperator r) {
        // Invoke operator= to assign *this operator
        RankTwoOperator<3, 3> &h = (*this);
        h[0][0] = (1.0 / 4.0) * (r[1](r[1]) + r[2](r[2]));
        h[0][1] = (-1.0 / 4.0) * r[0](r[1]);
        h[0][2] = (-1.0 / 4.0) * r[0](r[2]);
        h[1][0] = (-1.0 / 4.0) * r[1](r[0]);
        h[1][1] = (1.0 / 4.0) * (r[0](r[0]) + r[2](r[2]));
        h[1][2] = (-1.0 / 4.0) * r[1](r[2]);
        h[2][0] = (-1.0 / 4.0) * r[2](r[0]);
        h[2][1] = (-1.0 / 4.0) * r[2](r[1]);
        h[2][2] = (1.0 / 4.0) * (r[0](r[0]) + r[1](r[1]));
        h[0][0].name() = "h_BB_dia[x,x]";
        h[0][1].name() = "h_BB_dia[x,y]";
        h[0][2].name() = "h_BB_dia[x,z]";
        h[1][0].name() = "h_BB_dia[y,x]";
        h[1][1].name() = "h_BB_dia[y,y]";
        h[1][2].name() = "h_BB_dia[y,z]";
        h[2][0].name() = "h_BB_dia[z,x]";
        h[2][1].name() = "h_BB_dia[z,y]";
        h[2][2].name() = "h_BB_dia[z,z]";
    }
};

} // namespace mrchem
