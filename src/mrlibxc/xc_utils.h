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

#include <Eigen/Core>
#include <MRCPP/MWFunctions>
#include <XCFun/xcfun.h>

namespace mrdft {
namespace xc_utils {

Eigen::MatrixXi build_output_mask(bool is_lda, bool is_spin_sep, int order);
Eigen::VectorXi build_density_mask(bool is_lda, bool is_spin_sep, int order);

std::vector<mrcpp::FunctionNode<3> *> fetch_nodes(int n, mrcpp::FunctionTreeVector<3> &inp);
Eigen::MatrixXd compress_nodes(std::vector<mrcpp::FunctionNode<3> *> &inp_nodes);
void expand_nodes(std::vector<mrcpp::FunctionNode<3> *> &out_nodes, Eigen::MatrixXd &out_data);

mrcpp::FunctionTreeVector<3> log_gradient(mrcpp::DerivativeOperator<3> &diff_oper, mrcpp::FunctionTree<3> &rho);

} // namespace xc_utils
} // namespace mrdft
