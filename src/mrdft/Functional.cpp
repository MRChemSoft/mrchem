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
#include "Factory.h"

namespace mrdft {



void Functional::set_libxc_functional_object(std::vector<xc_func_type> libxc_objects_, std::vector<double> libxc_coeffs_) {
    libxc = Factory::libxc;
    libxc_objects = libxc_objects_;
    libxc_coeffs  = libxc_coeffs_;
    if (libxc) {
        std::cout << "RUNNING LIBXC" << std::endl;
    } else {
        std::cout << "RUNNING XCFUN" << std::endl;
    }
}


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
        // NB: the data is stored colomn major, i.e. two consecutive points of for example energy density, are not consecutive in memory
        // That means that we cannot extract the energy density data with out.row(0).data() for example.
        if (calc) xcfun_eval(xcfun.get(), inp.col(i).data(), out.col(i).data());
    }


    return out;
}


/** @brief Run a collection of grid points through XCFun
 *
 * Each column corresponds to one grid point.
 * From a performance point of view, (in pre and postprocessing) it is much more
 * efficient to have the two consecutive points in two consecutive adresses in memory
 *
 * param[in] inp_data Matrix of input values
 * param[out] out_data Matrix of output values
 */
Eigen::MatrixXd Functional::evaluate_transposed(Eigen::MatrixXd &inp) const {
    int nInp = xcfun_input_length(xcfun.get());  // Input parameters to XCFun
    int nOut = xcfun_output_length(xcfun.get()); // Input parameters to XCFun
    int nPts = inp.rows();
    if (nInp != inp.cols()) MSG_ABORT("Invalid input");
    // std::cout << "nInp: " << nInp << ", nOut: " << nOut << std::endl;
    
    Eigen::MatrixXd rho_spin  = Eigen::MatrixXd::Zero(nPts, 1);

    Eigen::MatrixXd out       = Eigen::MatrixXd::Zero(nPts, nOut);
    Eigen::MatrixXd out_libxc = Eigen::MatrixXd::Zero(nPts, nOut);
    Eigen::VectorXd out_sxc   = Eigen::VectorXd::Zero(nPts);
    Eigen::VectorXd exc, vxc, sxc, xxc, yxc, zxc, sigma, inp_row, out_row;



    if (Factory::libxc) {
        for (size_t i = 0; i < libxc_objects.size(); i++) {
            switch (libxc_objects[i].info->family) {
                case XC_FAMILY_LDA:
                case XC_FAMILY_HYB_LDA:
                        exc = Eigen::VectorXd::Zero(nPts);
                        vxc = Eigen::VectorXd::Zero(nPts);
                    if (isSpin()){
                        for (size_t k = 0; k < nPts; k++) {
                            rho_spin(k*2, 0)   = inp(k, 0);
                            rho_spin(k*2+1, 0) = inp(k, 1);
                        }
                        std::cout << "CONSTRUCTS SPIN RHO " << std::endl;
                        xc_lda_exc_vxc(&libxc_objects[i], nPts, rho_spin.col(0).data(), exc.data(), vxc.data());
                        for (size_t j = 0; j < nPts; ++j) {
                            //  xcfun computes rho * exc for energy density, so we do the same
                            //    aka xcfun calculates actual energy density while libxc calculates 
                            //    energy density per electron density
                            out_libxc(j, 0) += exc[j] * libxc_coeffs[i] * inp(j, 0);
                            out_libxc(j, 1) += vxc[j] * libxc_coeffs[i];
                        }
                    } else {
                        xc_lda_exc_vxc(&libxc_objects[i], nPts, inp.col(0).data(), exc.data(), vxc.data());
                        for (size_t j = 0; j < nPts; ++j) {
                            //  xcfun computes rho * exc for energy density, so we do the same
                            //    aka xcfun calculates actual energy density while libxc calculates 
                            //    energy density per electron density
                            out_libxc(j, 0) += exc[j] * libxc_coeffs[i] * inp(j, 0);
                            out_libxc(j, 1) += vxc[j] * libxc_coeffs[i];
                        }
                    }
                    break;
                case XC_FAMILY_GGA:
                case XC_FAMILY_HYB_GGA:
                    exc   = Eigen::VectorXd::Zero(nPts);
                    vxc   = Eigen::VectorXd::Zero(nPts);
                    sxc   = Eigen::VectorXd::Zero(nPts);
                    xxc   = Eigen::VectorXd::Zero(nPts);
                    yxc   = Eigen::VectorXd::Zero(nPts);
                    zxc   = Eigen::VectorXd::Zero(nPts);
                    sigma = Eigen::VectorXd::Zero(nPts);
                    sigma = inp.col(1) * inp.col(1) + inp.col(2) * inp.col(2) + inp.col(3) * inp.col(3);
                    xc_gga_exc_vxc(&libxc_objects[i], nPts, inp.col(0).data(), sigma.data(),
                        exc.data(), vxc.data(), sxc.data());
                    // Convert sxc to derivatives of x, y and z so output are similar to xcfun
                    xxc = 2 * sxc * inp.col(1);
                    yxc = 2 * sxc * inp.col(2);
                    zxc = 2 * sxc * inp.col(3);
                    for (size_t j = 0; j < nPts; ++j) {
                        //  xcfun computes rho * exc for energy density, so we do the same
                        //    aka xcfun calculates actual energy density while libxc calculates 
                        //    energy density per electron density
                        out_libxc(j, 0) += exc[j] * libxc_coeffs[i] * inp(j, 0);
                        out_libxc(j, 1) += vxc[j] * libxc_coeffs[i];
                        out_libxc(j, 2) += xxc[j] * libxc_coeffs[i];
                        out_libxc(j, 3) += yxc[j] * libxc_coeffs[i];
                        out_libxc(j, 4) += zxc[j] * libxc_coeffs[i];
                        out_sxc[j]      += sxc[j]; // only for debugging purposes
                    }
                    break;
                default:
                break;
            }
        }
    } else {
        inp_row = Eigen::VectorXd::Zero(nInp);
        out_row = Eigen::VectorXd::Zero(nOut);
        for (int i = 0; i < nPts; i++) {
            bool calc = true;
            if (isSpin()) {
                if (inp(i, 0) < cutoff and inp(i, 1) < cutoff) calc = false;
            } else {
                if (inp(i, 0) < cutoff) calc = false;
            }
            for (int j = 0; j < nInp; j++) inp_row(j) = inp(i, j);
            if (calc) xcfun_eval(xcfun.get(), inp_row.data(), out_row.data());
            for (int j = 0; j < nOut; j++) out(i, j) = out_row(j); 
            // if (calc) { // This to make out be 0 when calc = false?
            //     for (int j = 0; j < nInp; j++) inp_row(j) = inp(i, j);
            //     xcfun_eval(xcfun.get(), inp_row.data(), out_row.data());
            //     for (int j = 0; j < nOut; j++) out(i, j) = out_row(j);
            // }
        }
    }


    // for (size_t k = 0; k < nInp; k++) {
    //     for (size_t l = 0; l < nPts; l++) {
    //         if (abs(inp(l, k) - inp_row(k)) < 1e-10 && abs(inp_row(k)) > 1e-15) {
    //             std::cout << "XCFun inp row index: " << k << " = Libxc inp matrix index (row, col): " << l << ", " << k << " XCFun inp: " << inp_row(k) << std::endl;
    //         } 
    //     }
    // }

    // std::cout << std::endl;

    // std::cout << "----------------------" << std::endl << 
    // "XCFun inp row           : " << inp_row(0) << std::endl <<
    // "LibXC inp row inp(0, 0) : " << inp(0, 0) << std::endl <<
    // "LibXC inp row inp(0, 1) : " << inp(0, 1) << std::endl <<
    // "LibXC inp row inp(0, 2) : " << inp(0, 2) << std::endl <<
    // "LibXC inp row inp(1, 0) : " << inp(1, 0) << std::endl <<
    // "LibXC inp row inp(1, 1) : " << inp(1, 1) << std::endl <<
    // "LibXC inp row inp(1, 2) : " << inp(1, 2) << std::endl <<
    // "LibXC inp row inp(2, 0) : " << inp(2, 0) << std::endl <<
    // "LibXC inp row inp(2, 1) : " << inp(2, 1) << std::endl <<
    // "LibXC inp row inp(2, 2) : " << inp(2, 2) << std::endl <<
    // "----------------------" << std::endl;

    


    // // !! Debug !!

    // if (nInp == 1) {
    //     if (libxc) return out_libxc;
    //     if (!libxc) return out;
    // }

    // double dx = 0.0, dy = 0.0, dz = 0.0;

    // for (size_t i = 0; i < nPts; i++) {
    //     double sxc_x_xcfun = out(i, 2) / (2 * inp(i, 1));
    //     double sxc_y_xcfun = out(i, 3) / (2 * inp(i, 2));
    //     double sxc_z_xcfun = out(i, 4) / (2 * inp(i, 3));

    //     dx += sxc_x_xcfun - out_sxc[i];
    //     dy += sxc_y_xcfun - out_sxc[i];
    //     dz += sxc_z_xcfun - out_sxc[i];
    // }

    // dx /= nPts;
    // dy /= nPts;
    // dz /= nPts;

    // std::cout << "Cutoff: " << cutoff << std::endl;
    // std::cout << "SXC avg dev: " << dx << " " << dy << " " << dz << std::endl;

    // for (size_t i = 0; i < nOut; i++) {
    //     double max_dev = 0.0, mean_dev = 0.0;
    //     size_t max_i = 0;

    //     for (size_t j = 0; j < nPts; j++) {
    //         double dev = abs(out(j, i) - out_libxc(j, i));

    //         mean_dev += dev;

    //         if (dev > max_dev) {
    //             max_dev = dev;
    //             max_i = j;
    //         }
    //     }

    //     mean_dev /= nPts;

    //     std::cout << "Iteration nr.: " << i << " max_dev: " << max_dev << " mean_dev: " << mean_dev << std::endl;
    //     std::cout << "Current max i: " << max_i << std::endl <<
    //     "----------------------" << std::endl <<
    //     "Input: " << std::endl <<
    //     "    Density   : " << inp(max_i, 0) << std::endl <<
    //     "    dr/dx     : " << inp(max_i, 1) << std::endl <<
    //     "    dr/dy     : " << inp(max_i, 2) << std::endl <<
    //     "    dr/dz     : " << inp(max_i, 3) << std::endl <<
    //     "    sigma     : " << sigma[max_i] << std::endl <<
    //     "Output Libxc  : " << std::endl <<
    //     "    EXC       : " << out_libxc(max_i, 0) << std::endl <<
    //     "    VXC       : " << out_libxc(max_i, 1) << std::endl <<
    //     "    SXC       : " << sxc[max_i] << std::endl <<
    //     "    SXC (est) : " << out_libxc(max_i, 2) * out_libxc(max_i, 2) + out_libxc(max_i, 3) * out_libxc(max_i, 3) +
    //                           out_libxc(max_i, 4) * out_libxc(max_i, 4) << std::endl <<
    //     "    dEXC/dx   : " << out_libxc(max_i, 2) << std::endl <<
    //     "    dEXC/dy   : " << out_libxc(max_i, 3) << std::endl <<
    //     "    dEXC/dz   : " << out_libxc(max_i, 4) << std::endl <<
    //     "Output XCFun  : " << std::endl <<
    //     "    EXC       : " << out(max_i, 0) << std::endl <<
    //     "    VXC       : " << out(max_i, 1) << std::endl <<
    //     "    dEXC/dx   : " << out(max_i, 2) << std::endl <<
    //     "    dEXC/dy   : " << out(max_i, 3) << std::endl <<
    //     "    dEXC/dz   : " << out(max_i, 4) << std::endl <<
    //     "----------------------" << std::endl;
    // }
    // std::cout << std::endl;

    if (libxc) return out_libxc;
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

/** @brief Contract a collection of grid points
 *
 * Each column corresponds to one set of grid points.
 *
 * param[in] xc_data Matrix of functional partial derivative values
 * param[in] d_data Matrix of density input values
 * param[out] out_data Matrix of contracted output values
 */
Eigen::MatrixXd Functional::contract_transposed(Eigen::MatrixXd &xc_data, Eigen::MatrixXd &d_data) const {
    auto nPts = xc_data.rows();
    auto nFcs = getCtrOutputLength();
    Eigen::MatrixXd out_data = Eigen::MatrixXd::Zero(nPts, nFcs);
    out_data.col(0) = xc_data.col(0); // we always keep the energy functional

    for (int i = 0; i < this->xc_mask.rows(); i++) {
        Eigen::VectorXd cont_i = Eigen::VectorXd::Zero(nPts);
        for (int j = 0; j < this->xc_mask.cols(); j++) {
            Eigen::VectorXd cont_ij = Eigen::VectorXd::Zero(nPts);
            int xc_idx = this->xc_mask(i, j);
            int d_idx = this->d_mask(j);
            if (d_idx >= 0) {
                //elementwise product of one column of xc_data and d_data
                out_data.col(i + 1) += xc_data.col(xc_idx).cwiseProduct(d_data.col(d_idx));
            } else {
                out_data.col(i + 1) += xc_data.col(xc_idx);
            }
        }
    }
    return out_data;
}


/** @brief  Evaluates XC functional and derivatives for a given NodeIndex
 *
 * The electronic densities (total/alpha/beta) are given as input.
 * The values of the zero order densities and their gradient are sent to xcfun.
 * The output of xcfun must then be combined ("contract") with the gradients
 * of the higher order densities.
 *
 * XCFunctional output (with k=1 and explicit derivatives):
 *
 * LDA: \f$ \left(F_{xc}, \frac{\partial F_{xc}}{\partial \rho}\right) \f$
 *
 * GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho},
 *  \frac{\partial F_{xc}}{\partial \rho_x},
 *  \frac{\partial F_{xc}}{\partial \rho_y},
 *  \frac{\partial F_{xc}}{\partial \rho_z}\right) \f$
 *
 * Spin LDA: \f$ \left(F_{xc}, \frac{\partial F_{xc}}{\partial \rho^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho^\beta}\right) \f$
 *
 * Spin GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_x^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_y^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_z^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho_x^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_y^\beta},
 *  \frac{\partial F_{xc}}{\partial \rho_z^\beta}
 *  \right) \f$
 *
 * XCFunctional output (with k=1 and gamma-type derivatives):
 *
 * GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho},
 *  \frac{\partial F_{xc}}{\partial \gamma} \f$
 *
 * Spin GGA: \f$ \left(F_{xc},
 *  \frac{\partial F_{xc}}{\partial \rho^\alpha},
 *  \frac{\partial F_{xc}}{\partial \rho^\beta },
 *  \frac{\partial F_{xc}}{\partial \gamma^{\alpha \alpha}},
 *  \frac{\partial F_{xc}}{\partial \gamma^{\alpha \beta }},
 *  \frac{\partial F_{xc}}{\partial \gamma^{\beta  \beta }}
 *  \right) \f$
 *
 * param[in] inp Input values
 * param[out] xcNodes Output values
 *
 */
void Functional::makepot(mrcpp::FunctionTreeVector<3> &inp, std::vector<mrcpp::FunctionNode<3> *> xcNodes)  const {
    if (this->log_grad){
        MSG_ERROR("log_grad not implemented");
    }

    mrcpp::NodeIndex<3> nodeIdx = xcNodes[0]->getNodeIndex();
    mrcpp::FunctionTree<3>* rho0=std::get<1>(inp[0]);
    mrcpp::MWNode<3> node(rho0->getNode(nodeIdx),true,false); //copy node from rho, but do not copy coef
    int ncoefs = rho0->getTDim() * rho0->getKp1_d();
    int xcfun_inpsize = 1; // rho
    int spinsize = 1; // paired
    if (isSpin()) spinsize = 2; // alpha, beta
    xcfun_inpsize *= spinsize; // alpha and beta
    if (isGGA()) xcfun_inpsize *= 4; // add gradient (3 components for each spin)

    Eigen::MatrixXd xcfun_inp(ncoefs, xcfun_inpsize); //input for xcfun
    double* coef = node.getCoefs();

    for (int i = 0; i < spinsize; i++) {
        // make cv representation of density
        mrcpp::FunctionTree<3>* rho=std::get<1>(inp[i]);
        // we link into the node, in order to be able to do a mwtransform without copying the data back and forth
        node.attachCoefs(xcfun_inp.col(i).data());
        for (int j = 0; j < ncoefs; j++) xcfun_inp(j,i) = rho->getNode(nodeIdx).getCoefs()[j];
        node.mwTransform(mrcpp::Reconstruction);
        node.cvTransform(mrcpp::Forward);

        if (isGGA()) {
            //make gradient of input
            for (int d = 0; d < 3; d++) {
                node.attachCoefs(xcfun_inp.col(spinsize + 3*i + d).data());

                mrcpp::DerivativeCalculator<3> derivcalc(d, *this->derivOp, *rho);
                // derive rho and put result into xcfun_inp aka node
                derivcalc.calcNode(rho->getNode(nodeIdx), node);
                // make cv representation of gradient of density
                node.mwTransform(mrcpp::Reconstruction);
                node.cvTransform(mrcpp::Forward);
            }
       }
    }

    // send rho and grad rho to xcfun
    Eigen::MatrixXd xc_out = Functional::evaluate_transposed(xcfun_inp);

    // make gradient of the higher order densities
    //order:
    // rho_a_1
    // rho_b_1
    // drho_a_1/dx
    // drho_a_1/dy
    // drho_a_1/dz
    // drho_b_1/dx
    // drho_b_1/dy
    // drho_b_1/dz
    int ctrsize = inp.size()-spinsize; //number of higher order inputs
    int d_datasize = ctrsize;
    if (isGGA()) d_datasize *= 4; // add gradient (3 components for each higher order rho)
    Eigen::MatrixXd d_data = Eigen::MatrixXd::Zero(ncoefs, d_datasize);
    if (d_datasize > 0) {
        for (int i = 0; i < ctrsize; i++) {
            // make cv representation of density
            mrcpp::FunctionTree<3>* rho = std::get<1>(inp[i+spinsize]);
            // we link into the node, in order to be able to do a mwtransform without copying the data back and forth
            node.attachCoefs(d_data.col(i).data());
            for (int j = 0; j < ncoefs; j++) d_data(j,i) = rho->getNode(nodeIdx).getCoefs()[j];
            node.mwTransform(mrcpp::Reconstruction);
            node.cvTransform(mrcpp::Forward);
            if (isGGA()) {
                //make gradient of input
                for (int d = 0; d < 3; d++) {
                    node.attachCoefs(d_data.col(ctrsize + 3*i + d).data());
                    mrcpp::DerivativeCalculator<3> derivcalc(d, *this->derivOp, *rho);
                    derivcalc.calcNode(rho->getNode(nodeIdx), node);
                    // make cv representation of gradient of density
                    node.mwTransform(mrcpp::Reconstruction);
                    node.cvTransform(mrcpp::Forward);
                }
            }
        }
    }

    Eigen::MatrixXd Ctrout = contract_transposed(xc_out, d_data); //size output: LDA=1, GGA=4, spin *2

    // postprocess
    //For SpinGGA:
    //f_xc         : out[0] = inp[0]
    //df_xc/drho_a : out[1] = inp[1] - div(inp[3,4,5])
    //df_xc/drho_b : out[2] = inp[2] - div(inp[6,7,8])
    int xc_outsize = 2;
    if (isSpin()) xc_outsize = 3;
    for (int i = 0; i < xc_outsize; i++) {
        // from cv to node values
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
                //xcNodes[i] = Ctrout[i] - div(Ctrout[d_i])
                for (int j = 0; j < ncoefs; j++) xcNodes[i]->getCoefs()[j] -= noded.getCoefs()[j];
            }
        }
    }
    node.attachCoefs(coef); // restablish the original link (for proper destructor behaviour)
}
} // namespace mrdft