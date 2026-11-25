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

#include "Functional.h"

#include <MRCPP/Printer>

namespace mrdft {

void Functional::print_functional_references() const {
    if (mrcpp::mpi::world_rank != 0) return;
    int outfile_txt_width = 75;
    xclib->printFunctionalReference(outfile_txt_width);
    if (isMetaGGA()) println(0, " Meta-GGA lower tau cutoff (relative to tau_unif): " << lower_tau_cutoff);
}

double Functional::amountEXX() const {
    return xclib->getAmountExx();
}

void Functional::evaluate_data(const Eigen::MatrixXd &inp, Eigen::MatrixXd &out) const {
    // Should input be object instead of matrix so we can do inp.gradient()
    // and not deal with specific rows and columns?
    // would also fix problem with knowing length of input matrix
    // could also return functional data object
    int nInp = numIn();
    int nOut = numOut();
    size_t nPts = inp.cols();
    if (nInp != inp.rows()) {
        std::ostringstream oss;
        oss << "Invalid input: expected matrix with " << nInp << " rows, got " << inp.rows() << "!\n";
        MSG_ABORT(oss.str());
    }
    if (nOut != out.rows()) {
        std::ostringstream oss;
        oss << "Invalid output: expected matrix with " << nOut << " rows, got " << out.rows() << "!\n";
        MSG_ABORT(oss.str());
    }
    out.setZero();

    xclib->callLibEval(inp, out, nPts);
}

Eigen::MatrixXd Functional::evaluate(Eigen::MatrixXd &inp) const {
    // For efficiency: transpose inp and out matrices
    // NB: the data is stored column major, i.e. two consecutive points of for example energy density are not consecutive in memory.
    // That means that we cannot extract the energy density data with out.row(0).data() for example.
    Eigen::MatrixXd inp_trans(inp.transpose());
    Eigen::MatrixXd out_trans(numOut(), inp.rows());
    evaluate_data(inp_trans, out_trans);
    return out_trans.transpose();
}

Eigen::MatrixXd Functional::contract(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const {
    auto nPts = xc_data.rows();
    auto nFcs = getCtrOutputLength();
    Eigen::MatrixXd out_data = Eigen::MatrixXd::Zero(nPts, nFcs);
    out_data.col(0) = xc_data.col(0); // we always keep the energy functional

    for (int i = 0; i < this->xc_mask.rows(); i++) {
        for (int j = 0; j < this->xc_mask.cols(); j++) {
            int xc_idx = this->xc_mask(i, j);
            int d_idx = this->d_mask(j);
            if (d_idx >= 0) {
                // elementwise product of one column of xc_data and d_data
                out_data.col(i + 1) += xc_data.col(xc_idx).cwiseProduct(d_data.col(d_idx));
            } else {
                out_data.col(i + 1) += xc_data.col(xc_idx);
            }
        }
    }
    return out_data;
}

void Functional::makepot(mrcpp::FunctionTreeVector<3> &inp, std::vector<mrcpp::FunctionNode<3> *> xcNodes) const {
    if (this->log_grad) { MSG_ERROR("log_grad not implemented"); } // see Functional::log_grad

    mrcpp::NodeIndex<3> nodeIdx = xcNodes[0]->getNodeIndex();
    mrcpp::FunctionTree<3> *rho0 = std::get<1>(inp[0]);
    mrcpp::MWNode<3> node(rho0->getNode(nodeIdx), true, false); // copy node from rho, but do not copy coef
    int ncoefs = rho0->getTDim() * rho0->getKp1_d();
    int spinsize = densityChannels();
    int xclib_inpsize = numIn();

    Eigen::MatrixXd xclib_inp(ncoefs, xclib_inpsize); // input for the XC library
    // Meta-GGA (unpolarized): inp is (rho, tau) or, if XCPotentialD1 supplied grad rho built from the orbital
    // gradients, (rho, tau, d/dx rho, d/dy rho, d/dz rho). Then that grad rho is used instead of differentiating
    // rho, so that it is consistent with tau (built from the same orbital gradients).
    const bool orbital_grad_rho = isMetaGGA() && !isSpin() && static_cast<int>(inp.size()) >= spinsize + 4;
    double *coef = node.getCoefs();

    for (int i = 0; i < spinsize; i++) {
        // make cv representation of density
        mrcpp::FunctionTree<3> *rho = std::get<1>(inp[i]);
        auto &rhoNode = rho->getNode(nodeIdx);
        // we link into the node, in order to be able to do a mwtransform without copying the data back and forth
        node.attachCoefs(xclib_inp.col(i).data());
        for (int j = 0; j < ncoefs; j++) xclib_inp(j, i) = rhoNode.getCoefs()[j];
        node.mwTransform(mrcpp::Reconstruction);
        node.cvTransform(mrcpp::Forward);

        if (usesGradients()) {
            // make gradient of input
            for (int d = 0; d < 3; d++) {
                node.attachCoefs(xclib_inp.col(spinsize + 3 * i + d).data());
                if (orbital_grad_rho) {
                    mrcpp::FunctionTree<3> *grad_d = std::get<1>(inp[spinsize + 1 + d]);
                    for (int j = 0; j < ncoefs; j++) xclib_inp(j, spinsize + 3 * i + d) = grad_d->getNode(nodeIdx).getCoefs()[j];
                } else {
                    // derive rho and put result into xclib_inp aka node
                    mrcpp::DerivativeCalculator<3> derivcalc(d, *this->derivOp, *rho);
                    derivcalc.calcNode(rho->getNode(nodeIdx), node);
                }
                // make cv representation of gradient of density
                node.mwTransform(mrcpp::Reconstruction);
                node.cvTransform(mrcpp::Forward);
            }
        }
        if (isMetaGGA()) {
            int tau_col = spinsize * 4 + i;
            int tau_input_index = spinsize;

            if (tau_input_index >= static_cast<int>(inp.size())) {
                MSG_ABORT("Functional::makepot: tau input not found in FunctionTreeVector for meta-GGA.\n");
            }

            mrcpp::FunctionTree<3> *tau = std::get<1>(inp[tau_input_index]);
            if (!tau) {
                MSG_ABORT("Functional::makepot: tau pointer is null for meta-GGA.\n");
            }

            node.attachCoefs(xclib_inp.col(tau_col).data());
            for (int j = 0; j < ncoefs; j++) {
                xclib_inp(j, tau_col) = tau->getNode(nodeIdx).getCoefs()[j];
            }
            node.mwTransform(mrcpp::Reconstruction);
            node.cvTransform(mrcpp::Forward);
        }
    }

    // Meta-GGA input columns: 0 rho, 1-3 grad rho, 4 tau.
    // 1) Clamp tau >= tau_W (1 + 1e-10), tau_W = |grad rho|^2 / (8 rho): numerical noise can violate tau >= tau_W.
    //    The 1e-10 margin keeps Libxc's min(tau_W/tau, 1) on one derivative branch.
    // 2) Smooth lower cutoff (c = lower_tau_cutoff), tau_c = c tau_unif, tau_unif = 3/10 (3 pi^2)^(2/3) rho^(5/3):
    //    tau_eff = tau_c + tau^2 / (4 tau_c) where tau < 2 tau_c, otherwise tau. Keeps v_tau bounded where tau -> 0
    //    at finite rho (e.g. the H2 bond midpoint). The chain rule is applied to the potentials after evaluate().
    std::vector<double> smooth_x;  // x = tau/tau_c at points with x < 2, else -1 (not smoothed)
    std::vector<double> smooth_tc; // tau_c at those points
    if (isMetaGGA() && !isSpin()) {
        if (lower_tau_cutoff > 0.0) {
            smooth_x.assign(ncoefs, -1.0);
            smooth_tc.assign(ncoefs, 0.0);
        }
        for (int j = 0; j < ncoefs; j++) {
            double rho = xclib_inp(j, 0);
            if (rho <= 0.0) continue;
            double g2 = xclib_inp(j, 1) * xclib_inp(j, 1) + xclib_inp(j, 2) * xclib_inp(j, 2) + xclib_inp(j, 3) * xclib_inp(j, 3);
            double tau_min = g2 / (8.0 * rho) * (1.0 + 1.0e-10);
            if (xclib_inp(j, 4) < tau_min) xclib_inp(j, 4) = tau_min;
            if (lower_tau_cutoff > 0.0) {
                double tau_c = lower_tau_cutoff * 2.871234000188191 * std::pow(rho, 5.0 / 3.0);
                double tau = xclib_inp(j, 4);
                if (tau < 2.0 * tau_c) {
                    smooth_x[j] = tau / tau_c;
                    smooth_tc[j] = tau_c;
                    xclib_inp(j, 4) = tau_c + tau * tau / (4.0 * tau_c);
                }
            }
        }
    }

    // send rho, grad rho (and tau) to the XC library
    Eigen::MatrixXd xc_out = Functional::evaluate(xclib_inp);

    // Chain rule for the smooth cutoff (unpolarized meta-GGA output columns: 0 f, 1 dF/drho,
    // 2-4 dF/d(grad rho), 5 dF/dtau): with tau_eff(rho, tau),
    //   v_tau = dF/dtau_eff * x/2,   v_rho += dF/dtau_eff * (1 - x^2/4) * (5/3) tau_c / rho.
    if (!smooth_x.empty()) {
        for (int j = 0; j < ncoefs; j++) {
            double x = smooth_x[j];
            if (x < 0.0) continue;
            double rho = xclib_inp(j, 0);
            double f_tau = xc_out(j, 5);
            xc_out(j, 1) += f_tau * (1.0 - 0.25 * x * x) * (5.0 / 3.0) * smooth_tc[j] / rho;
            xc_out(j, 5) = f_tau * 0.5 * x;
        }
    }

    // Meta-GGA: screen the low-density tail, where the functional is evaluated on numerical
    // noise (rho < 1e-15 gave v_tau ~ 1e10-1e24 and SCF divergence). Rows are grid points.
    if (isMetaGGA() && !isSpin()) {
        const double rho_min = std::max(this->cutoff, 1.0e-10);
        for (int j = 0; j < ncoefs; j++) {
            if (xclib_inp(j, 0) < rho_min) xc_out.row(j).setZero();
        }
    }

    // make gradient of the higher order densities
    // layout per higher-order rho when gradients are used:
    //   1 density value + 3 spatial gradient components (x, y, z)
    // order in inp / d_data:
    //   rho_a_1
    //   rho_b_1
    //   drho_a_1/dx
    //   drho_a_1/dy
    //   drho_a_1/dz
    //   drho_b_1/dx
    //   drho_b_1/dy
    //   drho_b_1/dz
    // Meta-GGA: the inputs after rho are tau and grad rho, not higher-order densities
    // (meta-GGA response is not implemented, XCPotentialD2 aborts).
    int ctrsize = isMetaGGA() ? 0 : inp.size() - spinsize; // number of higher-order inputs (each: density + 3 gradient components if usesGradients())
    int d_datasize = ctrsize;
    if (usesGradients()) d_datasize *= 4; // 1 density + 3 gradient components for each higher-order rho
    Eigen::MatrixXd d_data = Eigen::MatrixXd::Zero(ncoefs, d_datasize);
    if (d_datasize > 0) {
        for (int i = 0; i < ctrsize; i++) {
            // make cv representation of density
            mrcpp::FunctionTree<3> *rho = std::get<1>(inp[i + spinsize]);
            // we link into the node, in order to be able to do a mwtransform without copying the data back and forth
            node.attachCoefs(d_data.col(i).data());
            for (int j = 0; j < ncoefs; j++) d_data(j, i) = rho->getNode(nodeIdx).getCoefs()[j];
            node.mwTransform(mrcpp::Reconstruction);
            node.cvTransform(mrcpp::Forward);
            if (usesGradients()) {
                // make gradient of input
                for (int d = 0; d < 3; d++) {
                    node.attachCoefs(d_data.col(ctrsize + 3 * i + d).data());
                    mrcpp::DerivativeCalculator<3> derivcalc(d, *this->derivOp, *rho);
                    derivcalc.calcNode(rho->getNode(nodeIdx), node);
                    // make cv representation of gradient of density
                    node.mwTransform(mrcpp::Reconstruction);
                    node.cvTransform(mrcpp::Forward);
                }
            }
        }
    }

    // Meta-GGA (order 1 only): the library output is already (f, dF/drho, dF/d(grad rho) (3), dF/dtau);
    // there is nothing to contract
    Eigen::MatrixXd Ctrout = isMetaGGA() ? xc_out : contract(xc_out, d_data); // size output: LDA=1, GGA=4, spin *2

    // postprocess
    // For spin-polarized GGA:
    // f_xc         : out[0] = inp[0]
    // df_xc/drho_a : out[1] = inp[1] - div(inp[3,4,5])
    // df_xc/drho_b : out[2] = inp[2] - div(inp[6,7,8])
    // Meta-GGA (unpolarized) has the GGA layout in columns 0-4, plus v_tau = dF/dtau in column 5 (-> xcNodes[2])
    int xc_outsize = densityChannels() + 1;
    for (int i = 0; i < xc_outsize; i++) {
        // from cv to node values
        node.attachCoefs(Ctrout.col(i).data());
        node.cvTransform(mrcpp::Backward);
        node.mwTransform(mrcpp::Compression);
        for (int j = 0; j < ncoefs; j++) xcNodes[i]->getCoefs()[j] = Ctrout(j, i);
        xcNodes[i]->setHasCoefs();
        if (usesGradients() and i > 0) {
            for (int d = 0; d < 3; d++) {
                node.attachCoefs(Ctrout.col(xc_outsize + 3 * (i - 1) + d).data());
                node.cvTransform(mrcpp::Backward);
                node.mwTransform(mrcpp::Compression);
                node.calcNorms();
                mrcpp::DerivativeCalculator<3> derivcalc(d, *this->derivOp, *rho0); // TODO: define outside loops
                mrcpp::MWNode<3> noded(rho0->getNode(nodeIdx), true, false);
                derivcalc.calcNode(node, noded);
                // xcNodes[i] = Ctrout[i] - div(Ctrout[d_i])
                for (int j = 0; j < ncoefs; j++) xcNodes[i]->getCoefs()[j] -= noded.getCoefs()[j];
            }
        }
    }
    if (isMetaGGA()) {
        if (xcNodes.size() < 3) MSG_ABORT("Functional::makepot: meta-GGA output vector too small for tau potential.");
        node.attachCoefs(Ctrout.col(5).data());
        node.cvTransform(mrcpp::Backward);
        node.mwTransform(mrcpp::Compression);
        for (int j = 0; j < ncoefs; j++) xcNodes[2]->getCoefs()[j] = Ctrout(j, 5);
        xcNodes[2]->setHasCoefs();
    }
    node.attachCoefs(coef); // reestablish the original link (for proper destructor behaviour)
}
} // namespace mrdft
