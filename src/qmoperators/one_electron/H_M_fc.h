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

#include "chemistry/PhysicalConstants.h"

#include "DeltaOperator.h"
#include "H_B_spin.h"

namespace mrchem {

/** @class H_M_fc
 *
 * @brief Fermi-Contact operator
 *
 * Interaction operator obtained by differentiating the spin Hamiltonian wrt
 * the nuclear magnetic moment of nucleus K:
 *
 * dH/dM_K = H_M_pso + H_M_sd + H_M_fc
 *
 * H_M_fc = -\frac{8\pi\alpha^2}{3} \sum_j \delta(r_{jK})m_j
 *
 * where m_j is the magnetic moment of the electron.
 */

class H_M_fc final : public RankOneOperator<3> {
public:
    H_M_fc(const mrcpp::Coord<3> &o, double proj_prec, double smooth_prec = -1.0)
            : H_M_fc(H_B_spin(), DeltaOperator(o, proj_prec, smooth_prec)) {}

    H_M_fc(H_B_spin s, DeltaOperator delta) {
        const double coef = -(8.0 / 3.0) * mrcpp::pi;
        const double alpha_2 = PhysicalConstants::get("fine_structure_constant") * PhysicalConstants::get("fine_structure_constant") * 1000000.0;

        // Invoke operator= to assign *this operator
        RankOneOperator<3> &h = (*this);
        h[0] = (coef * alpha_2) * delta * s[0];
        h[1] = (coef * alpha_2) * delta * s[1];
        h[2] = (coef * alpha_2) * delta * s[2];
        h[0].name() = "h_M_fc[x]";
        h[1].name() = "h_M_fc[y]";
        h[2].name() = "h_M_fc[z]";
    }
};

} // namespace mrchem
