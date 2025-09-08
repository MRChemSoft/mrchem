/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) ...
 */

#include <MRCPP/Printer>

#include "Functional.h"

namespace mrdft {

namespace {
// Put XCFun into a valid “user eval” mode so input_length/output_length work.
inline void ensure_xcfun_user_setup(xcfun_t* xf, int order, bool spin, bool is_gga) {
    const unsigned int mode      = 1u;                    // partial derivatives
    const unsigned int func_type = is_gga ? 1u : 0u;      // 0=LDA, 1=GGA
    const unsigned int dens_type = spin ? 2u : 1u;        // 1 (unpol) or 2 (pol)
    const unsigned int laplacian = 0u;
    const unsigned int kinetic   = 0u;
    const unsigned int current   = 0u;
    const unsigned int exp_deriv = is_gga ? 1u : 0u;      // explicit derivs for GGA only

    xcfun_user_eval_setup(xf, order, func_type, dens_type,
                          mode, laplacian, kinetic, current, exp_deriv);
}
} // namespace

/** @brief Run a collection of grid points through XCFun
 *
 * Each row corresponds to one grid point.
 */
Eigen::MatrixXd Functional::evaluate(Eigen::MatrixXd &inp) const {
    // Make sure XCFun knows what shape to expect/produce
    ensure_xcfun_user_setup(xcfun.get(), order, isSpin(), isGGA());

    int nInp = xcfun_input_length(xcfun.get());   // input parameters to XCFun
    int nOut = xcfun_output_length(xcfun.get());  // output parameters from XCFun
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

/** @brief Run a collection of grid points through XCFun
 *
 * Each column corresponds to one grid point.
 */
Eigen::MatrixXd Functional::evaluate_transposed(Eigen::MatrixXd &inp) const {
    // Make sure XCFun knows what shape to expect/produce
    ensure_xcfun_user_setup(xcfun.get(), order, isSpin(), isGGA());

    int nInp = xcfun_input_length(xcfun.get());   // input parameters to XCFun
    int nOut = xcfun_output_length(xcfun.get());  // output parameters from XCFun
    int nPts = inp.rows();
    if (nInp != inp.cols()) MSG_ABORT("Invalid input");

    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(nPts, nOut);
    Eigen::VectorXd inp_row = Eigen::VectorXd::Zero(nInp);
    Eigen::VectorXd out_row = Eigen::VectorXd::Zero(nOut);
    for (int i = 0; i < nPts; i++) {
        bool calc = true;
        if (isSpin()) {
            if (inp(i, 0) < cutoff and inp(i, 1) < cutoff) calc = false;
        } else {
            if (inp(i, 0) < cutoff) calc = false;
        }
        if (calc) {
            inp_row = inp.row(i);
            xcfun_eval(xcfun.get(), inp_row.data(), out_row.data());
            out.row(i) = out_row;
        }
    }
    return out;
}

/** @brief Contract (row-major) */
Eigen::MatrixXd Functional::contract(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const {
    auto nPts = xc_data.cols();
    auto nFcs = getCtrOutputLength();
    Eigen::MatrixXd out_data = Eigen::MatrixXd::Zero(nFcs, nPts);
    out_data.row(0) = xc_data.row(0); // keep energy functional

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
        out_data.row(i + 1) = cont_i;
    }
    return out_data;
}

/** @brief Contract (column-major) */
Eigen::MatrixXd Functional::contract_transposed(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const {
    auto nPts = xc_data.rows();
    auto nFcs = getCtrOutputLength();
    Eigen::MatrixXd out_data = Eigen::MatrixXd::Zero(nPts, nFcs);
    out_data.col(0) = xc_data.col(0); // keep energy functional

    for (int i = 0; i < this->xc_mask.rows(); i++) {
        for (int j = 0; j < this->xc_mask.cols(); j++) {
            int xc_idx = this->xc_mask(i, j);
            int d_idx = this->d_mask(j);
            if (d_idx >= 0) {
                out_data.col(i + 1) += xc_data.col(xc_idx).cwiseProduct(d_data.col(d_idx));
            } else {
                out_data.col(i + 1) += xc_data.col(xc_idx);
            }
        }
    }
    return out_data;
}

/** @brief makepot: compute XC and contract */
void Functional::makepot(mrcpp::FunctionTreeVector<3> &inp, std::vector<mrcpp::FunctionNode<3> *> xcNodes)  const {
    if (this->log_grad){
        MSG_ERROR("log_grad not implemented");
    }

    mrcpp::NodeIndex<3> nodeIdx = xcNodes[0]->getNodeIndex();
    mrcpp::FunctionTree<3>* rho0=std::get<1>(inp[0]);
    mrcpp::MWNode<3> node(rho0->getNode(nodeIdx),true,false);
    int ncoefs = rho0->getTDim() * rho0->getKp1_d();
    int xcfun_inpsize = 1; // rho
    int spinsize = 1; // paired
    if (isSpin()) spinsize = 2; // alpha, beta
    xcfun_inpsize *= spinsize; // alpha and beta
    if (isGGA()) xcfun_inpsize *= 4; // add gradient (3 components for each spin)

    Eigen::MatrixXd xcfun_inp(ncoefs, xcfun_inpsize);
    double* coef = node.getCoefs();

    for (int i = 0; i < spinsize; i++) {
        mrcpp::FunctionTree<3>* rho=std::get<1>(inp[i]);
        node.attachCoefs(xcfun_inp.col(i).data());
        for (int j = 0; j < ncoefs; j++) xcfun_inp(j,i) = rho->getNode(nodeIdx).getCoefs()[j];
        node.mwTransform(mrcpp::Reconstruction);
        node.cvTransform(mrcpp::Forward);

        if (isGGA()) {
            for (int d = 0; d < 3; d++) {
                node.attachCoefs(xcfun_inp.col(spinsize + 3*i + d).data());
                mrcpp::DerivativeCalculator<3> derivcalc(d, *this->derivOp, *rho);
                derivcalc.calcNode(rho->getNode(nodeIdx), node);
                node.mwTransform(mrcpp::Reconstruction);
                node.cvTransform(mrcpp::Forward);
            }
       }
    }

    // NB: VIRTUAL DISPATCH here (so LibXC adapters can override)!
    Eigen::MatrixXd xc_out = this->evaluate_transposed(xcfun_inp);

    int ctrsize = inp.size()-spinsize;
    int d_datasize = ctrsize;
    if (isGGA()) d_datasize *= 4;
    Eigen::MatrixXd d_data = Eigen::MatrixXd::Zero(ncoefs, d_datasize);
    if (d_datasize > 0) {
        for (int i = 0; i < ctrsize; i++) {
            mrcpp::FunctionTree<3>* rho = std::get<1>(inp[i+spinsize]);
            node.attachCoefs(d_data.col(i).data());
            for (int j = 0; j < ncoefs; j++) d_data(j,i) = rho->getNode(nodeIdx).getCoefs()[j];
            node.mwTransform(mrcpp::Reconstruction);
            node.cvTransform(mrcpp::Forward);
            if (isGGA()) {
                for (int d = 0; d < 3; d++) {
                    node.attachCoefs(d_data.col(ctrsize + 3*i + d).data());
                    mrcpp::DerivativeCalculator<3> derivcalc(d, *this->derivOp, *rho);
                    derivcalc.calcNode(rho->getNode(nodeIdx), node);
                    node.mwTransform(mrcpp::Reconstruction);
                    node.cvTransform(mrcpp::Forward);
                }
            }
        }
    }

    Eigen::MatrixXd Ctrout = contract_transposed(xc_out, d_data);

    int xc_outsize = 2;
    if (isSpin()) xc_outsize = 3;
    for (int i = 0; i < xc_outsize; i++) {
        node.attachCoefs(Ctrout.col(i).data());
        node.cvTransform(mrcpp::Backward);
        node.mwTransform(mrcpp::Compression);
        for (int j = 0; j < ncoefs; j++) xcNodes[i]->getCoefs()[j] = Ctrout(j,i);
        xcNodes[i]->setHasCoefs();
        if (isGGA() and i>0) {
            for (int d = 0; d < 3; d++) {
                node.attachCoefs(Ctrout.col(xc_outsize + 3*(i-1) + d).data());
                node.cvTransform(mrcpp::Backward);
                node.mwTransform(mrcpp::Compression);
                node.calcNorms();
                mrcpp::DerivativeCalculator<3> derivcalc(d,*this->derivOp, *rho0);//TODO: define outside loops
                mrcpp::MWNode<3> noded(rho0->getNode(nodeIdx),true,false);
                derivcalc.calcNode(node, noded);
                for (int j = 0; j < ncoefs; j++) xcNodes[i]->getCoefs()[j] -= noded.getCoefs()[j];
            }
        }
    }
    node.attachCoefs(coef);
}

} // namespace mrdft
