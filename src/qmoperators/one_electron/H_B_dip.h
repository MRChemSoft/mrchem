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

#include "AngularMomentumOperator.h"

namespace mrchem {

/** @class H_B_dip
 *
 * @brief Magnetic dipole operator
 *
 * Interaction operator obtained by differentiating the spin Hamiltonian wrt
 * the external magnetic field B:
 *
 * dH/dB = H_B_dip + H_B_spin
 *
 * H_B_dip = \frac{1}{2} \sum_j l_{jO}
 *
 * where l_{jO} is the orbital angular momentum.
 */

class H_B_dip final : public RankOneOperator<3> {
public:
    H_B_dip(std::shared_ptr<mrcpp::DerivativeOperator<3>> D, const mrcpp::Coord<3> &o)
            : H_B_dip(AngularMomentumOperator(D, o)) {}

    explicit H_B_dip(AngularMomentumOperator l) {
        // Invoke operator= to assign *this operator
        RankOneOperator<3> &h = (*this);
        h[0] = 0.5 * l[0];
        h[1] = 0.5 * l[1];
        h[2] = 0.5 * l[2];
        h[0].name() = "h_B_dip[x]";
        h[1].name() = "h_B_dip[y]";
        h[2].name() = "h_B_dip[z]";
    }
};

} // namespace mrchem
