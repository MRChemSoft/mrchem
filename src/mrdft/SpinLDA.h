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
#include "MRCPP/Printer"
#include "xc_utils.h"

namespace mrdft {

/**
 * @class LDA
 * @brief Local Density Approximation xc functional
 * @details Implements a spin‑polarized LDA functional on the MRA grid
 */
class SpinLDA final : public Functional {
public:
    /**
     * @brief Construct a LDA functional
     * @param[in] k Polynomial order of the functional derivatives (0: energy,
     *             1: potential, etc.)
     * @param[in] f XC library handle (Libxc/XCFun); ownership is transferred
     */
    SpinLDA(int k, XCLib_p &f)
            : Functional(k, f) {
        xc_mask = xc_utils::build_output_mask(true, true, this->order);
        d_mask = xc_utils::build_density_mask(true, true, this->order);
    }
    ~SpinLDA() override = default;

    bool isSpin() const override { return true; }
    bool isGGA() const override { return false; }
    bool isMetaGGA() const override { return false; }

    int numIn() const override { return 2; }                 ///< @brief Number of input components: 1 densities (alpha, beta)
    int numOut() const override { return xclib->getnOut(); } ///< @brief Number of raw outputs provided by the xc backend
    int densityChannels() const override { return 2; }       ///< @brief Two density channels; rho_alpha, rho_beta
    bool usesGradients() const override { return false; }    ///< @brief LDA does not require density gradients

private:
    mrcpp::FunctionTreeVector<3> rho_a;
    mrcpp::FunctionTreeVector<3> rho_b;

    int getCtrOutputLength() const override { return 3; } ///< @brief Number of contracted outputs (energy + 2 components for LDA)

    /** @brief Clear internal functions
     *
     * Ownership of densities (alpha, beta) is outside MRDFT -> clear
     */
    void clear() override {
        mrcpp::clear(this->rho_a, false);
        mrcpp::clear(this->rho_b, false);
    }
};

} // namespace mrdft
