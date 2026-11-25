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
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <vector>

#include "Factory.h"
#include "utils/print_utils.h"

namespace mrdft {

void Functional::print_functional_references() const {
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
    auto nPts = xc_data.cols();
    auto nFcs = getCtrOutputLength();
    Eigen::MatrixXd out_data = Eigen::MatrixXd::Zero(nFcs, nPts);
    out_data.row(0) = xc_data.row(0); // we always keep the energy functional

    for (int i = 0; i < this->xc_mask.rows(); i++) {
        // Eigen::VectorXd cont_i = Eigen::VectorXd::Zero(nPts);
        for (int j = 0; j < this->xc_mask.cols(); j++) {
            int xc_idx = this->xc_mask(i, j);
            int d_idx = this->d_mask(j);
            if (d_idx >= 0) {
                out_data.row(i + 1).array() += xc_data.row(xc_idx).array() * d_data.row(d_idx).array();
            } else {
                out_data.row(i + 1).array() += xc_data.row(xc_idx).array();
            }
        }
        // out_data.row(i + 1) = cont_i; // The first column contains the energy functional
    }
    return out_data;
}

Eigen::MatrixXd Functional::contract_transposed(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const {
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
    if (this->log_grad) { MSG_ERROR("log_grad not implemented"); }

    mrcpp::NodeIndex<3> nodeIdx = xcNodes[0]->getNodeIndex();
    mrcpp::FunctionTree<3> *rho0 = std::get<1>(inp[0]);
    mrcpp::MWNode<3> node(rho0->getNode(nodeIdx), true, false); // copy node from rho, but do not copy coef
    int ncoefs = rho0->getTDim() * rho0->getKp1_d();
    int spinsize = densityChannels();
    int xclib_inpsize = spinsize * (usesGradients() ? 4 : 1);
    if (isMetaGGA()) xclib_inpsize += spinsize;

    // int xcfun_inpsize = 1; // rho
    // int spinsize = isSpin() ? 2 : 1; // 1 for unpolarized, 2 for polarized
    // xcfun_inpsize *= spinsize; // alpha and beta
    // if (isGGA() || isMetaGGA()) xcfun_inpsize *= 4; // add gradient (3 components for each spin)
    // if(isMetaGGA()) xcfun_inpsize += spinsize;      // add tau (1 for each spin)

    Eigen::MatrixXd xclib_inp(ncoefs, xclib_inpsize); // input for xcfun
    // Meta-GGA (unpolarized): if XCPotentialD1 supplied grad rho built from the orbital gradients
    // (inputs 2, 3, 4 after rho and tau), use it instead of differentiating rho, so that it is
    // consistent with tau (built from the same orbital gradients).
    const bool orbital_grad_rho = isMetaGGA() && !isSpin() && static_cast<int>(inp.size()) >= spinsize + 4;
    double *coef = node.getCoefs();

    for (int i = 0; i < spinsize; i++) {
        // make cv representation of density
        mrcpp::FunctionTree<3> *rho = std::get<1>(inp[i]);
        auto &rhoNode = rho->getNode(nodeIdx);
        node.attachCoefs(xclib_inp.col(i).data());
        for (int j = 0; j < ncoefs; j++) xclib_inp(j, i) = rhoNode.getCoefs()[j];
        node.mwTransform(mrcpp::Reconstruction);
        node.cvTransform(mrcpp::Forward);

        if (usesGradients()) {
            for (int d = 0; d < 3; d++) {
                node.attachCoefs(xclib_inp.col(spinsize + 3 * i + d).data());
                if (orbital_grad_rho) {
                    mrcpp::FunctionTree<3> *grad_d = std::get<1>(inp[spinsize + 1 + d]);
                    for (int j = 0; j < ncoefs; j++) xclib_inp(j, spinsize + 3 * i + d) = grad_d->getNode(nodeIdx).getCoefs()[j];
                } else {
                    mrcpp::DerivativeCalculator<3> derivcalc(d, *this->derivOp, *rho);
                    derivcalc.calcNode(rho->getNode(nodeIdx), node);
                }
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

    int ctrsize = inp.size() - spinsize;
    int d_datasize = ctrsize;
    if (usesGradients()) d_datasize *= 4;
    if (isMetaGGA() && !isSpin()) {
        ctrsize = 0;
        d_datasize = 0;
    } else if (isMetaGGA() && isSpin()) {
        MSG_ABORT("Functional::makepot: spin-polarized meta-GGA not supported.");
    }
    if (ctrsize != getCtrInputLength()) {
        std::ostringstream oss;
        oss << "Functional::makepot: mismatch between ctr inputs ("
            << ctrsize << ") and expected (" << getCtrInputLength() << ").\n" << " spinsize: " << spinsize << " inp.size(): " << inp.size();
        MSG_ABORT(oss.str());
    }

    Eigen::MatrixXd d_data = Eigen::MatrixXd::Zero(ncoefs, d_datasize);
    if (d_datasize > 0) {
        for (int i = 0; i < ctrsize; i++) {
            mrcpp::FunctionTree<3> *rho = std::get<1>(inp[i + spinsize]);
            node.attachCoefs(d_data.col(i).data());
            for (int j = 0; j < ncoefs; j++) d_data(j, i) = rho->getNode(nodeIdx).getCoefs()[j];
            node.mwTransform(mrcpp::Reconstruction);
            node.cvTransform(mrcpp::Forward);
            if (usesGradients()) {
                for (int d = 0; d < 3; d++) {
                    node.attachCoefs(d_data.col(ctrsize + 3 * i + d).data());
                    mrcpp::DerivativeCalculator<3> derivcalc(d, *this->derivOp, *rho);
                    derivcalc.calcNode(rho->getNode(nodeIdx), node);
                    node.mwTransform(mrcpp::Reconstruction);
                    node.cvTransform(mrcpp::Forward);
                }
            }
        }
    }

    Eigen::MatrixXd Ctrout;
    if (isMetaGGA() && !isSpin() && getCtrInputLength() == 0) {
        Ctrout = xc_out;
    } else {
        Ctrout = contract_transposed(xc_out, d_data);
    }

    if (isMetaGGA() && !isSpin() && getCtrInputLength() == 0) {
        if (xcNodes.size() < 3) {
            MSG_ABORT("Functional::makepot: meta-GGA output vector too small for tau potential.\n");
        }

        node.attachCoefs(Ctrout.col(0).data());
        node.cvTransform(mrcpp::Backward);
        node.mwTransform(mrcpp::Compression);
        for (int j = 0; j < ncoefs; j++) xcNodes[0]->getCoefs()[j] = Ctrout(j, 0);
        xcNodes[0]->setHasCoefs();

        node.attachCoefs(Ctrout.col(1).data());
        node.cvTransform(mrcpp::Backward);
        node.mwTransform(mrcpp::Compression);
        for (int j = 0; j < ncoefs; j++) xcNodes[1]->getCoefs()[j] = Ctrout(j, 1);
        xcNodes[1]->setHasCoefs();
        // Gradient part of the potential, as in the GGA branch below:
        // v_xc = dF/drho - div( dF/d(grad rho) ), with dF/d(d_d rho) in columns 2, 3, 4.
        for (int d = 0; d < 3; d++) {
            node.attachCoefs(Ctrout.col(2 + d).data());
            node.cvTransform(mrcpp::Backward);
            node.mwTransform(mrcpp::Compression);
            node.calcNorms();
            mrcpp::DerivativeCalculator<3> derivcalc(d, *this->derivOp, *rho0);
            mrcpp::MWNode<3> noded(rho0->getNode(nodeIdx), true, false);
            derivcalc.calcNode(node, noded);
            for (int j = 0; j < ncoefs; j++) xcNodes[1]->getCoefs()[j] -= noded.getCoefs()[j];
        }

        node.attachCoefs(Ctrout.col(5).data());
        node.cvTransform(mrcpp::Backward);
        node.mwTransform(mrcpp::Compression);
        for (int j = 0; j < ncoefs; j++) xcNodes[2]->getCoefs()[j] = Ctrout(j, 5);
        xcNodes[2]->setHasCoefs();
    } else {
        int xc_outsize = densityChannels() + 1;
        for (int i = 0; i < xc_outsize; i++) {
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
                    mrcpp::DerivativeCalculator<3> derivcalc(d, *this->derivOp, *rho0);
                    mrcpp::MWNode<3> noded(rho0->getNode(nodeIdx), true, false);
                    derivcalc.calcNode(node, noded);
                    for (int j = 0; j < ncoefs; j++) xcNodes[i]->getCoefs()[j] -= noded.getCoefs()[j];
                }
            }
        }
    }
    node.attachCoefs(coef);
}
} // namespace mrdft
