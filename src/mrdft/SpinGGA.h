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
 * @class spinGGA
 * @brief Generalized Gradient Approximation xc functional
 * @details Implements a spin‑polarized GGA functional on the MRA grid
 */
class SpinGGA final : public Functional {
public:
    /**
     * @brief Construct a spin GGA functional
     * @param[in] k Polynomial order of the functional derivatives (0: energy,
     *             1: potential, etc.)
     * @param[in] f XC library handle (Libxc/XCFun); ownership is transferred
     * @param[in] d Numerical derivative operator used to compute density gradients
     */
    SpinGGA(int k, XCLib_p &f, std::shared_ptr<mrcpp::DerivativeOperator<3>> &d)
            : Functional(k, f)
            , derivative(std::move(d)) {
        xc_mask = xc_utils::build_output_mask(false, true, this->order);
        d_mask = xc_utils::build_density_mask(false, true, this->order);
    }
    ~SpinGGA() override = default;

    bool isSpin() const override { return true; }
    bool isGGA() const override { return true; }
    bool isMetaGGA() const override { return false; }

    int numIn() const override { return 8; }                 ///< @brief Number of input components: 2 densities + 6 gradient components
    int numOut() const override { return xclib->getnOut(); } ///< @brief Number of raw outputs provided by the xc backend
    int densityChannels() const override { return 2; }       ///< @brief Two density channels; rho_alpha, rho_beta
    bool usesGradients() const override { return true; }     ///< @brief spinGGA requires density gradients

private:
    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative{nullptr};
    mrcpp::FunctionTreeVector<3> rho_a;
    mrcpp::FunctionTreeVector<3> rho_b;
    mrcpp::FunctionTreeVector<3> grad_a;
    mrcpp::FunctionTreeVector<3> grad_b;

    int getCtrOutputLength() const override { return 9; } ///< @brief Number of contracted outputs (energy + 8 components for spinGGA)

    /** @brief Clear internal functions
     *
     * Ownership of densities (alpha, beta) is outside MRDFT -> clear
     * Ownership of gradients (alpha, beta) is inside MRDFT -> free
     */
    void clear() override {
        mrcpp::clear(this->rho_a, false);
        mrcpp::clear(this->rho_b, false);
        mrcpp::clear(this->grad_a, true);
        mrcpp::clear(this->grad_b, true);
    }
};

} // namespace mrdft
