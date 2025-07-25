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

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/trees/FunctionNode.h>

#include "xc_utils.h"

namespace mrdft {

namespace xc_utils {
void fill_output_mask(Eigen::MatrixXi &mask, int start);
} // namespace xc_utils

Eigen::MatrixXi xc_utils::build_output_mask(bool is_lda, bool is_spin_sep, int order) {
    int start = 2;
    bool is_gga = not is_lda;
    Eigen::MatrixXi mask(1, 1);
    mask << 1;
    switch (order) {
        case 0:
            break;
        case 1:
            if (is_lda and is_spin_sep) {
                mask.resize(2, 1);
                mask << 1, 2;
            } else if (is_gga and not is_spin_sep) {
                mask.resize(4, 1);
                mask << 1, 2, 3, 4;
            } else if (is_gga and is_spin_sep) {
                mask.resize(8, 1);
                mask << 1, 2, 3, 4, 5, 6, 7, 8;
            }
            break;
        case 2:
            if (is_lda and is_spin_sep) {
                start = 3;
                mask.resize(2, 2);
            } else if (is_gga and not is_spin_sep) {
                start = 5;
                mask.resize(4, 4);
            } else if (is_gga and is_spin_sep) {
                start = 9;
                mask.resize(8, 8);
            }
            fill_output_mask(mask, start);
            break;
        default:
            MSG_ABORT("Not implemented");
    }
    return mask;
}

Eigen::VectorXi xc_utils::build_density_mask(bool is_lda, bool is_spin_sep, int order) {
    bool is_gga = not is_lda;
    Eigen::VectorXi mask(1);
    switch (order) {
        case 0:
        case 1:
            mask(0) = -1;
            break;
        case 2:
            mask(0) = 0;
            if (is_lda and is_spin_sep) {
                mask.resize(2);
                mask << 0, 1;
            } else if (is_gga and not is_spin_sep) {
                mask.resize(4);
                mask << 0, 1, 2, 3;
            } else if (is_gga and is_spin_sep) {
                mask.resize(8);
                mask << 0, 1, 2, 3, 4, 5, 6, 7;
            }
            break;
        default:
            MSG_ABORT("Not implemented");
    }
    return mask;
}

void xc_utils::fill_output_mask(Eigen::MatrixXi &mask, int value) {
    for (int i = 0; i < mask.rows(); i++) {
        mask(i, i) = value;
        value++;
        for (int j = i + 1; j < mask.cols(); j++) {
            mask(i, j) = value;
            mask(j, i) = value;
            value++;
        }
    }
}

/** @brief Fetch specific node from several FunctionTrees
 *
 * This will retrieve one node from each of the input trees and put them
 * into a vector of FunctionNodes. The node is fetched from position n
 * in the respective endNodeTables, which means that the tree structures
 * must be identical for this routine to work as intended.
 *
 * param[in] n Node position in EndNodeTable
 * param[in] inp_trees Array of FunctionTrees
 * param[out] out_nodes Array of FunctionNodes
 */
std::vector<mrcpp::FunctionNode<3> *> xc_utils::fetch_nodes(int n, mrcpp::FunctionTreeVector<3> &inp_trees) {
    std::vector<mrcpp::FunctionNode<3> *> out_nodes;
    for (auto i = 0; i < inp_trees.size(); i++) {
        auto &iTree = mrcpp::get_func(inp_trees, i);
        auto &iNode = iTree.getEndFuncNode(n);
        out_nodes.push_back(&iNode);
    }
    return out_nodes;
}

/** @brief Collect data from FunctionNodes into a matrix
 *
 * Collects function values from the input nodes into the rows
 * of a matrix. Matrix dimension: rows = nNodes, cols = nCoefs.
 *
 * param[in] inp_nodes Array of FunctionNodes
 * param[out] out_data Matrix of function values
 */
Eigen::MatrixXd xc_utils::compress_nodes(std::vector<mrcpp::FunctionNode<3> *> &inp_nodes) {
    Eigen::MatrixXd out_data;
    auto nNodes = inp_nodes.size();
    if (nNodes > 0) {
        auto nCoefs = inp_nodes[0]->getNCoefs();
        out_data = Eigen::MatrixXd::Zero(nNodes, nCoefs);
        for (auto i = 0; i < nNodes; i++) {
            auto &node = inp_nodes[i];
            Eigen::VectorXd row_i;
            node->getValues(row_i);
            if (row_i.size() != nCoefs) MSG_ABORT("Size mismatch");
            out_data.row(i) = row_i;
        }
    }
    return out_data;
}

/** @brief Put data from a matrix into FunctionNodes
 *
 * Each row of the input data is used as function values
 * of the corresponding FunctionNode in the output vector.
 * Matrix dimension: rows = nNodes, cols = nCoefs.
 *
 * param[inout] out_nodes Array of FunctionNodes
 * param[in] inp_data Matrix of function values
 */
void xc_utils::expand_nodes(std::vector<mrcpp::FunctionNode<3> *> &out_nodes, Eigen::MatrixXd &inp_data) {
    auto nFuncs = out_nodes.size();
    if (inp_data.rows() != nFuncs) MSG_ERROR("Size mismatch " << inp_data.rows() << " vs " << nFuncs);

    for (auto i = 0; i < nFuncs; i++) {
        auto &node = out_nodes[i];
        node->setValues(inp_data.row(i));
    }
}

/** @brief Compute the gradient using a log parametrization
 *
 * zeta = log(inp_func)
 * grad(inp_func) = inp_func * grad(zeta)
 *
 * param[in] diff_oper Derivative operator
 * param[in] inp_func Function to differentiate
 * param[out] out_grad Gradient of input function
 */
mrcpp::FunctionTreeVector<3> xc_utils::log_gradient(mrcpp::DerivativeOperator<3> &diff_oper, mrcpp::FunctionTree<3> &inp_func) {
    mrcpp::FunctionTree<3> zeta(inp_func.getMRA());
    mrcpp::copy_grid(zeta, inp_func);
    mrcpp::copy_func(zeta, inp_func);
    for (auto i = 0; i < zeta.getNEndNodes(); i++) {
        auto &node = zeta.getEndFuncNode(i);
        Eigen::VectorXd values;
        node.getValues(values);
        for (auto j = 0; j < node.getNCoefs(); j++) {
            if (values[j] > mrcpp::MachineZero) {
                values[j] = std::log(values[j]);
            } else {
                values[j] = mrcpp::MachineZero;
            }
        }
        node.setValues(values);
    }
    zeta.mwTransform(mrcpp::BottomUp);

    mrcpp::FunctionTreeVector<3> grad_zeta = mrcpp::gradient(diff_oper, zeta);

    mrcpp::FunctionTreeVector<3> out_grad;
    for (int i = 0; i < 3; i++) {
        mrcpp::FunctionTree<3> *grad_comp = new mrcpp::FunctionTree<3>(inp_func.getMRA());
        mrcpp::copy_grid(*grad_comp, inp_func);
        mrcpp::multiply(-1.0, *grad_comp, 1.0, inp_func, mrcpp::get_func(grad_zeta, i));
        out_grad.push_back(std::make_tuple(1.0, grad_comp));
    }
    mrcpp::clear(grad_zeta, true);
    return out_grad;
}

} // namespace mrdft
