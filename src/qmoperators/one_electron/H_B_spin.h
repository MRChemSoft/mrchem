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

#include "chemistry/PhysicalConstants.h"
#include "tensor/RankOneOperator.h"

#include "SpinOperator.h"

namespace mrchem {

/** @class H_B_spin
 *
 * @brief Magnetic spin operator
 *
 * Interaction operator obtained by differentiating the spin Hamiltonian wrt
 * the external magnetic field B:
 *
 * dH/dB = H_B_dip + H_B_spin
 *
 * H_B_spin = -\sum_j m_j
 *
 * where m_j is the magnetic moment of the electron.
 */

class H_B_spin final : public RankOneOperator<3> {
public:
    H_B_spin()
            : H_B_spin(SpinOperator()) {}

    explicit H_B_spin(SpinOperator s) {
        const double g_e = PhysicalConstants::get("electron_g_factor");

        // Invoke operator= to assign *this operator
        RankOneOperator<3> &h = (*this);
        h[0] = (g_e / 2.0) * s[0];
        h[1] = (g_e / 2.0) * s[1];
        h[2] = (g_e / 2.0) * s[2];
        h[0].name() = "h_B_spin[x]";
        h[1].name() = "h_B_spin[y]";
        h[2].name() = "h_B_spin[z]";
    }
};

} // namespace mrchem
