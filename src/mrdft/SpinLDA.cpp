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

#include "MRCPP/MWFunctions"
#include "MRCPP/Printer"

#include "SpinLDA.h"
#include "xc_utils.h"

namespace mrdft {

SpinLDA::SpinLDA(int k, XC_p &f)
        : Functional(k, f) {
    xc_mask = xc_utils::build_output_mask(true, true, this->order);
    d_mask = xc_utils::build_density_mask(true, true, this->order);
}

/** @brief Clear internal functions
 *
 * Ownership of densities is outside MRDFT -> clear
 * Ownership of gradients is inside MRDFT -> free
 */
void SpinLDA::clear() {
    mrcpp::clear(this->rho_a, false);
    mrcpp::clear(this->rho_b, false);
}

/** @brief Number of function involved in contraction step */
int SpinLDA::getCtrInputLength() const {
    int length = -1;
    if (this->order < 2) length = 0;
    if (this->order == 2) length = 2;
    if (this->order > 2) NOT_IMPLEMENTED_ABORT;
    return length;
}

/** @brief Collect input functions to xcfun evaluation step
 *
 * For SpinLDA : [alpha_0, beta_0]
 */
mrcpp::FunctionTreeVector<3> SpinLDA::setupXCInput() {
    if (this->rho_a.size() < 1) MSG_ERROR("Alpha density not initialized");
    if (this->rho_b.size() < 1) MSG_ERROR("Beta density not initialized");

    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(this->rho_a[0]);
    out_vec.push_back(this->rho_b[0]);
    return out_vec;
}

/** @brief Collect input functions to contraction step
 *
 * For SpinLDA:
 * Ground State: No contraction, empty vector
 * Linear Response: [alpha_1, beta_1]
 * Higher Response: NOT_IMPLEMENTED
 */
mrcpp::FunctionTreeVector<3> SpinLDA::setupCtrInput() {
    if (this->order > 2) NOT_IMPLEMENTED_ABORT;
    mrcpp::FunctionTreeVector<3> out_vec;
    if (order == 2) {
        out_vec.push_back(this->rho_a[1]);
        out_vec.push_back(this->rho_b[1]);
    }
    return out_vec;
}

/** @brief Prepare input functions to xcfun
 *
 * Collects input densities and computes necessary gradients.
 *
 * Ordering of input:
 * inp_vec[0] = alpha_0
 * inp_vec[1] = beta_0
 * inp_vec[2] = alpha_1
 * inp_vec[3] = beta_1
 * ...
 */
void SpinLDA::preprocess(mrcpp::FunctionTreeVector<3> &inp_vec) {
    if (inp_vec.size() != 2 * this->order) MSG_ERROR("Invalid input length");
    if (this->rho_a.size() > 0) MSG_ERROR("Alpha density not empty");
    if (this->rho_b.size() > 0) MSG_ERROR("Beta density not empty");

    int n = 0;
    for (int i = 0; i < this->order; i++) {
        this->rho_a.push_back(inp_vec[n++]);
        this->rho_b.push_back(inp_vec[n++]);
    }
}

/** @brief Compute final output functions
 *
 * Combine the raw partial derivatives from xcfun into functional derivatives.
 *
 * For SpinLDA:
 * f_xc         : out[0] = inp[0]
 * df_xc/drho_a : out[1] = inp[1]
 * df_xc/drho_b : out[2] = inp[2]
 */
mrcpp::FunctionTreeVector<3> SpinLDA::postprocess(mrcpp::FunctionTreeVector<3> &inp_vec) {
    // Energy density
    mrcpp::FunctionTree<3> &f_xc = mrcpp::get_func(inp_vec, 0);
    inp_vec[0] = std::make_tuple<double, mrcpp::FunctionTree<3> *>(1.0, nullptr);

    // Alpha potential
    mrcpp::FunctionTree<3> &v_a = mrcpp::get_func(inp_vec, 1);
    inp_vec[1] = std::make_tuple<double, mrcpp::FunctionTree<3> *>(1.0, nullptr);

    // Beta potential
    mrcpp::FunctionTree<3> &v_b = mrcpp::get_func(inp_vec, 2);
    inp_vec[2] = std::make_tuple<double, mrcpp::FunctionTree<3> *>(1.0, nullptr);

    // Collect output
    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(std::make_tuple(1.0, &f_xc));
    out_vec.push_back(std::make_tuple(1.0, &v_a));
    out_vec.push_back(std::make_tuple(1.0, &v_b));
    return out_vec;
}
} // namespace mrdft
