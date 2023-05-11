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

#include <MRCPP/Printer>

#include "Functional.h"

namespace mrdft {

/** @brief Run a collection of grid points through XCFun
 *
 * Each row corresponds to one grid point.
 *
 * param[in] inp_data Matrix of input values
 * param[out] out_data Matrix of output values
 */
Eigen::MatrixXd Functional::evaluate(Eigen::MatrixXd &inp) const {
    int nInp = xcfun_input_length(xcfun.get());  // Input parameters to XCFun
    int nOut = xcfun_output_length(xcfun.get()); // Input parameters to XCFun
    int nPts = inp.cols();
    if (nInp != inp.rows()) MSG_ABORT("Invalid input");

    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(nOut, nPts);
    for (int i = 0; i < nPts; i++) {
        bool calc = true;
        if (isSpin()) {
            if (inp(0, i) < cutoff and inp(1, i) < cutoff) calc = false;
        } else {
            if (inp(0, i) < cutoff) calc = false;
        }
        if (calc) xcfun_eval(xcfun.get(), inp.col(i).data(), out.col(i).data());
    }
    return out;
}

/** @brief Contract a collection of grid points
 *
 * Each row corresponds to one grid point.
 *
 * param[in] xc_data Matrix of functional partial derivative values
 * param[in] d_data Matrix of density input values
 * param[out] out_data Matrix of contracted output values
 */
Eigen::MatrixXd Functional::contract(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const {
    auto nPts = xc_data.cols();
    auto nFcs = getCtrOutputLength();
    Eigen::MatrixXd out_data = Eigen::MatrixXd::Zero(nFcs, nPts);
    out_data.row(0) = xc_data.row(0); // we always keep the energy functional

    for (int i = 0; i < this->xc_mask.rows(); i++) {
        Eigen::VectorXd cont_i = Eigen::VectorXd::Zero(nPts);
        for (int j = 0; j < this->xc_mask.cols(); j++) {
            Eigen::VectorXd cont_ij = Eigen::VectorXd::Zero(nPts);
            int xc_idx = this->xc_mask(i, j);
            int d_idx = this->d_mask(j);
            if (d_idx >= 0) {
                cont_ij = xc_data.row(xc_idx).array() * d_data.row(d_idx).array();
            } else {
                cont_ij = xc_data.row(xc_idx);
            }
            cont_i += cont_ij;
        }
        out_data.row(i + 1) = cont_i; // The first column contains the energy functional
    }
    return out_data;
}

} // namespace mrdft
