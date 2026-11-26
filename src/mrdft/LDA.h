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

#include "Functional.h"
#include "xc_utils.h"

namespace mrdft {

/**
 * @class LDA
 * @brief Local Density Approximation xc functional
 * @details Implements the spin-restricted and the spin-polarized LDA functional on the MRA grid
 */
class LDA final : public Functional {
public:
    /**
     * @brief Construct a LDA functional
     * @param[in] k Polynomial order of the functional derivatives (0: energy,
     *             1: potential, etc.)
     * @param[in] s True for the spin-polarized case, false for spin-restricted
     * @param[in] f XC library handle (Libxc/XCFun); ownership is transferred
     */
    LDA(int k, bool s, XCLib_p &f)
            : Functional(k, s, f) {
        xc_mask = xc_utils::build_output_mask(true, spin, this->order);
        d_mask = xc_utils::build_density_mask(true, spin, this->order);
    }
    ~LDA() override = default;

    bool isGGA() const override { return false; }
    bool isMetaGGA() const override { return false; }

    int numOut() const override { return xclib->getnOut(); } ///< @brief Number of raw outputs provided by the xc backend
    bool usesGradients() const override { return false; }    ///< @brief LDA does not require density gradients
};

} // namespace mrdft
