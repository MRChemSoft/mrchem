/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "chemistry/PhysicalConstants.h"
#include "tensor/RankTwoOperator.h"

#include "NuclearGradientOperator.h"
#include "PositionOperator.h"

namespace mrchem {

/** @class H_BM_dia
 *
 * @brief Diamagnetic NMR operator
 *
 * Interaction operator obtained by differentiating the spin Hamiltonian wrt the
 * external magnetic field B and the nuclear magnetic moment M_K of nucleus K:
 *
 * d^2H/dBdM_K = H_BM_dia
 *
 * H_BM_dia = \frac{alpha^2}{2}
 *            \sum_j \frac{(r_{jO} \cdot r_{jK})1 - r_{jK}r_{jO}^T}{r_jK}^3}
 *
 * This operator is the transpose of H_MB_dia.
 */

class H_BM_dia final : public RankTwoOperator<3, 3> {
public:
    H_BM_dia(const mrcpp::Coord<3> &o, const mrcpp::Coord<3> &k, double proj_prec, double smooth_prec = -1.0)
            : H_BM_dia(PositionOperator(o), NuclearGradientOperator(1.0, k, proj_prec, smooth_prec)) {}

    H_BM_dia(PositionOperator r_o, NuclearGradientOperator r_rm3) {
        const double alpha_2 = PhysicalConstants::get("fine_structure_constant") * PhysicalConstants::get("fine_structure_constant") * 1000000.0;

        RankZeroOperator &o_x = r_o[0];
        RankZeroOperator &o_y = r_o[1];
        RankZeroOperator &o_z = r_o[2];
        RankZeroOperator &k_x = r_rm3[0];
        RankZeroOperator &k_y = r_rm3[1];
        RankZeroOperator &k_z = r_rm3[2];

        // Invoke operator= to assign *this operator
        RankTwoOperator<3, 3> &h = (*this);
        h[0][0] = (alpha_2 / 2.0) * (o_y(k_y) + o_z(k_z));
        h[0][1] = -(alpha_2 / 2.0) * o_y(k_x);
        h[0][2] = -(alpha_2 / 2.0) * o_z(k_x);
        h[1][0] = -(alpha_2 / 2.0) * o_x(k_y);
        h[1][1] = (alpha_2 / 2.0) * (o_x(k_x) + o_z(k_z));
        h[1][2] = -(alpha_2 / 2.0) * o_z(k_y);
        h[2][0] = -(alpha_2 / 2.0) * o_x(k_z);
        h[2][1] = -(alpha_2 / 2.0) * o_y(k_z);
        h[2][2] = (alpha_2 / 2.0) * (o_x(k_x) + o_y(k_y));
        h[0][0].name() = "h_BM_dia[x,x]";
        h[0][1].name() = "h_BM_dia[x,y]";
        h[0][2].name() = "h_BM_dia[x,z]";
        h[1][0].name() = "h_BM_dia[y,x]";
        h[1][1].name() = "h_BM_dia[y,y]";
        h[1][2].name() = "h_BM_dia[y,z]";
        h[2][0].name() = "h_BM_dia[z,x]";
        h[2][1].name() = "h_BM_dia[z,y]";
        h[2][2].name() = "h_BM_dia[z,z]";
    }
};

} // namespace mrchem
