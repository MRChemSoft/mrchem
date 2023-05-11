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

#include "MRCPP/MWOperators"
#include "MRCPP/Printer"

#include "SpinGGA.h"
#include "xc_utils.h"

namespace mrdft {

SpinGGA::SpinGGA(int k, XC_p &f, std::unique_ptr<mrcpp::DerivativeOperator<3>> &d)
        : Functional(k, f)
        , derivative(std::move(d)) {
    xc_mask = xc_utils::build_output_mask(false, true, this->order);
    d_mask = xc_utils::build_density_mask(false, true, this->order);
}

/** @brief Clear internal functions
 *
 * Ownership of densities is outside MRDFT -> clear
 * Ownership of gradients is inside MRDFT -> free
 */
void SpinGGA::clear() {
    mrcpp::clear(this->rho_a, false);
    mrcpp::clear(this->rho_b, false);
    mrcpp::clear(this->grad_a, true);
    mrcpp::clear(this->grad_b, true);
}

/** @brief Number of function involved in contraction step */
int SpinGGA::getCtrInputLength() const {
    int length = -1;
    if (this->order < 2) length = 0;
    if (this->order == 2) length = 8;
    if (this->order > 2) NOT_IMPLEMENTED_ABORT;
    return length;
}

/** @brief Collect input functions to xcfun evaluation step
 *
 * For SpinGGA : [alpha_0, beta_0, grad(alpha_0), grad(beta_0)]
 */
mrcpp::FunctionTreeVector<3> SpinGGA::setupXCInput() {
    if (this->rho_a.size() < 1) MSG_ERROR("Alpha density not initialized");
    if (this->rho_b.size() < 1) MSG_ERROR("Beta density not initialized");
    if (this->grad_a.size() < 3) MSG_ERROR("Alpha gradient not initialized");
    if (this->grad_b.size() < 3) MSG_ERROR("Beta gradient not initialized");

    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(this->rho_a[0]);
    out_vec.push_back(this->rho_b[0]);
    out_vec.insert(out_vec.end(), this->grad_a.begin(), this->grad_a.begin() + 3);
    out_vec.insert(out_vec.end(), this->grad_b.begin(), this->grad_b.begin() + 3);
    return out_vec;
}

/** @brief Collect input functions to contraction step
 *
 * For SpinGGA:
 * Ground State: No contraction, empty vector
 * Linear Response: [alpha_1, beta_1, grad(alpha_1), grad(beta_1)]
 * Higher Response: NOT_IMPLEMENTED
 */
mrcpp::FunctionTreeVector<3> SpinGGA::setupCtrInput() {
    if (this->order > 2) NOT_IMPLEMENTED_ABORT;
    mrcpp::FunctionTreeVector<3> out_vec;
    if (this->order == 2) {
        out_vec.push_back(this->rho_a[1]);
        out_vec.push_back(this->rho_b[1]);
        out_vec.insert(out_vec.end(), this->grad_a.begin() + 3, this->grad_a.begin() + 6);
        out_vec.insert(out_vec.end(), this->grad_b.begin() + 3, this->grad_b.begin() + 6);
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
void SpinGGA::preprocess(mrcpp::FunctionTreeVector<3> &inp_vec) {
    if (inp_vec.size() != 2 * this->order) MSG_ERROR("Invalid input length");
    if (this->rho_a.size() > 0) MSG_ERROR("Alpha density not empty");
    if (this->rho_b.size() > 0) MSG_ERROR("Beta density not empty");
    if (this->grad_a.size() > 0) MSG_ERROR("Alpha gradient not empty");
    if (this->grad_b.size() > 0) MSG_ERROR("Beta gradient not empty");

    int n = 0;
    for (int i = 0; i < this->order; i++) {
        this->rho_a.push_back(inp_vec[n++]);
        this->rho_b.push_back(inp_vec[n++]);
    }

    for (int i = 0; i < this->order; i++) {
        mrcpp::FunctionTreeVector<3> tmp_a, tmp_b;
        if (this->log_grad and i == 0) {
            tmp_a = xc_utils::log_gradient(*this->derivative, mrcpp::get_func(this->rho_a, i));
            tmp_b = xc_utils::log_gradient(*this->derivative, mrcpp::get_func(this->rho_b, i));
        } else {
            tmp_a = mrcpp::gradient(*this->derivative, mrcpp::get_func(this->rho_a, i));
            tmp_b = mrcpp::gradient(*this->derivative, mrcpp::get_func(this->rho_b, i));
        }
        this->grad_a.insert(this->grad_a.end(), tmp_a.begin(), tmp_a.end());
        this->grad_b.insert(this->grad_b.end(), tmp_b.begin(), tmp_b.end());
    }
}

/** @brief Compute final output functions
 *
 * Combine the raw partial derivatives from xcfun into functional derivatives.
 *
 * For SpinGGA:
 * f_xc         : out[0] = inp[0]
 * df_xc/drho_a : out[1] = inp[1] - div(inp[3,4,5])
 * df_xc/drho_b : out[2] = inp[2] - div(inp[6,7,8])
 */
mrcpp::FunctionTreeVector<3> SpinGGA::postprocess(mrcpp::FunctionTreeVector<3> &inp_vec) {
    // Energy density
    mrcpp::FunctionTree<3> &f_xc = mrcpp::get_func(inp_vec, 0);
    inp_vec[0] = std::make_tuple<double, mrcpp::FunctionTree<3> *>(1.0, nullptr);

    // Alpha part
    mrcpp::FunctionTree<3> &df_da = mrcpp::get_func(inp_vec, 1);
    mrcpp::FunctionTreeVector<3> df_dga(inp_vec.begin() + 3, inp_vec.begin() + 6);

    auto *tmp_a = new mrcpp::FunctionTree<3>(df_da.getMRA());
    mrcpp::divergence(*tmp_a, *this->derivative, df_dga);

    auto *v_a = new mrcpp::FunctionTree<3>(df_da.getMRA());
    mrcpp::build_grid(*v_a, df_da);
    mrcpp::build_grid(*v_a, *tmp_a);
    mrcpp::add(-1.0, *v_a, 1.0, df_da, -1.0, *tmp_a);
    delete tmp_a;

    // Beta part
    mrcpp::FunctionTree<3> &df_db = mrcpp::get_func(inp_vec, 2);
    mrcpp::FunctionTreeVector<3> df_dgb(inp_vec.begin() + 6, inp_vec.begin() + 9);

    auto *tmp_b = new mrcpp::FunctionTree<3>(df_db.getMRA());
    mrcpp::divergence(*tmp_b, *this->derivative, df_dgb);

    auto *v_b = new mrcpp::FunctionTree<3>(df_db.getMRA());
    mrcpp::build_grid(*v_b, df_db);
    mrcpp::build_grid(*v_b, *tmp_b);
    mrcpp::add(-1.0, *v_b, 1.0, df_db, -1.0, *tmp_b);
    delete tmp_b;

    // Collect output
    mrcpp::FunctionTreeVector<3> out_vec;
    out_vec.push_back(std::make_tuple(1.0, &f_xc));
    out_vec.push_back(std::make_tuple(1.0, v_a));
    out_vec.push_back(std::make_tuple(1.0, v_b));
    v_a = nullptr;
    v_b = nullptr;

    return out_vec;
}

} // namespace mrdft
