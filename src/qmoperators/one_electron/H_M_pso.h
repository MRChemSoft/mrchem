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
#include "tensor/RankOneOperator.h"

#include "MomentumOperator.h"
#include "NuclearGradientOperator.h"

namespace mrchem {

/** @class H_M_pso
 *
 * @brief Paramagnetic Spin-Orbit operator
 *
 * Interaction operator obtained by differentiating the spin Hamiltonian wrt
 * the nuclear magnetic moment of nucleus K:
 *
 * dH/dM_K = H_M_pso + H_M_sd + H_M_fc
 *
 * H_M_pso = alpha^2 \sum_j \frac{l_{jK}}{r_{jK}^3}
 *
 * where l_{jK} is the orbital angular momentum.
 */

class H_M_pso final : public RankOneOperator<3> {
public:
    H_M_pso(std::shared_ptr<mrcpp::DerivativeOperator<3>> D, const mrcpp::Coord<3> &k, double proj_prec, double smooth_prec = -1.0)
            : H_M_pso(MomentumOperator(D), NuclearGradientOperator(1.0, k, proj_prec, smooth_prec)) {}

    H_M_pso(MomentumOperator p, NuclearGradientOperator r_rm3) {
        const double alpha_2 = PhysicalConstants::get("fine_structure_constant") * PhysicalConstants::get("fine_structure_constant") * 1000000.0;

        RankZeroOperator &p_x = p[0];
        RankZeroOperator &p_y = p[1];
        RankZeroOperator &p_z = p[2];
        RankZeroOperator &x_rm3 = r_rm3[0];
        RankZeroOperator &y_rm3 = r_rm3[1];
        RankZeroOperator &z_rm3 = r_rm3[2];

        // Invoke operator= to assign *this operator
        RankOneOperator<3> &h = (*this);
        h[0] = alpha_2 * (y_rm3 * p_z - z_rm3 * p_y);
        h[1] = alpha_2 * (z_rm3 * p_x - x_rm3 * p_z);
        h[2] = alpha_2 * (x_rm3 * p_y - y_rm3 * p_x);
        h[0].name() = "h_M_pso[x]";
        h[1].name() = "h_M_pso[y]";
        h[2].name() = "h_M_pso[z]";
    }
};

} // namespace mrchem
