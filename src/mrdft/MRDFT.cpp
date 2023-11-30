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
#include "MRCPP/Timer"
#include <MRCPP/trees/FunctionNode.h>
#include <MRCPP/Parallel>

#include "Functional.h"
#include "MRDFT.h"
#include "xc_utils.h"

using namespace std;

namespace mrdft {

/** @brief Compute XC potentials from densities
 *
 * This routine computes the XC energy density and potentials on the
 * union grid of all the density input functions. Functional evaluation
 * and subsequent contraction is done node by node, to avoid explicit
 * construction of the huge number of intermediate functions.
 *
 * Ordering without spin:
 * inp_vec[0] = rho_0 (unperturbed)
 * inp_vec[1] = rho_1 (first order perturbed)
 * ...
 * out_vec[0] = f_xc (XC energy density)
 * out_vec[1] = v_xc (XC potential)
 *
 * Ordering with spin:
 * inp_vec[0] = alpha_0 (unperturbed)
 * inp_vec[1] = beta_0
 * inp_vec[2] = alpha_1 (first order perturbed)
 * inp_vec[3] = beta_1
 * ...
 * out_vec[0] = f_xc (XC energy density)
 * out_vec[1] = v_xc_a (XC alpha potential)
 * out_vec[2] = v_xc_b (XC beta potential)
 */
mrcpp::FunctionTreeVector<3> MRDFT::evaluate(mrcpp::FunctionTreeVector<3> &inp) {
    mrcpp::Timer t_tot, t_pre;
    grid().unify(inp);
    int nCoefs = mrcpp::get_func(inp, 0).getEndFuncNode(0).getNCoefs();
    int nOutCtr = functional().getCtrOutputLength();
    int nFcs = functional().getXCOutputLength();

    //Note: the other option will be removed after log_grad is implemented
    if (not functional().log_grad) {
        //Note: we do not need to make the complete energy density (PotVec[0]) for each mpi,
        //      instead we compute the energy directly.
        mrcpp::Timer t_post;
        int potvecSize = 2;
        if (functional().isSpin()) potvecSize = 3;
        mrcpp::FunctionTreeVector<3> PotVec = grid().generate(potvecSize);
        int nNodes = grid().size();
        int n_start = (mrcpp::mpi::wrk_rank * nNodes) / mrcpp::mpi::wrk_size;
        int n_end = ((mrcpp::mpi::wrk_rank + 1) * nNodes) / mrcpp::mpi::wrk_size;
        DoubleVector XCenergy = DoubleVector::Zero(1);
#pragma omp parallel
        {
#pragma omp for schedule(guided)
            for (int n = n_start; n < n_end; n++) {
               vector<mrcpp::FunctionNode<3> *> xcNodes = xc_utils::fetch_nodes(n, PotVec);
               functional().makepot(inp, xcNodes);
               XCenergy[0] += xcNodes[0]->integrate();
            }
        }
        // each mpi only has part of the results. All send their results to bank and then fetch
        if(mrcpp::mpi::wrk_size > 1) {
            // sum up the energy contrbutions from all mpi
            mrcpp::mpi::allreduce_vector(XCenergy, mrcpp::mpi::comm_wrk);
            // send to bank
            // note that this cannot be done in the omp loop above, because omp threads cannot use mpi
            mrcpp::BankAccount PotVecBank; // to put the PotVec
            for (int n = n_start; n < n_end; n++) {
                vector<mrcpp::FunctionNode<3> *> xcNodes = xc_utils::fetch_nodes(n, PotVec);
                //NB: we do not distribute the energy density (i=0). It is not used, since we have XCenergy
                for (int i = 1; i < potvecSize; i++)
                PotVecBank.put_data(potvecSize*n+i, nCoefs, xcNodes[i]->getCoefs());
            }
            for (int n = 0; n < nNodes; n++) {
                if(n >= n_start and n < n_end) continue; //no need to fetch own results
                vector<mrcpp::FunctionNode<3> *> xcNodes = xc_utils::fetch_nodes(n, PotVec);
                for (int i = 1; i < potvecSize; i++)
                PotVecBank.get_data(potvecSize*n+i, nCoefs, xcNodes[i]->getCoefs());
            }
        }
        functional().clear();
        int outNodes = 0;
        int outSize = 0;
        for (int i = 1; i < PotVec.size(); i++) {
            mrcpp::FunctionTree<3> &f_i = mrcpp::get_func(PotVec, i);
            f_i.mwTransform(mrcpp::BottomUp);
            f_i.calcSquareNorm();
            if(i>0)outNodes += f_i.getNNodes();
            if(i>0)outSize += f_i.getSizeNodes();
        }
        this->functional().XCenergy = XCenergy[0];
        mrcpp::print::tree(3, "Make potential", outNodes, outSize, t_post.elapsed());
        return PotVec;
    } else {

    functional().preprocess(inp);
    mrcpp::FunctionTreeVector<3> xcInpVec = functional().setupXCInput();
    mrcpp::FunctionTreeVector<3> ctrInpVec = functional().setupCtrInput();

    int inpNodes = 0;
    int inpSize = 0;
    for (auto i = 0; i < xcInpVec.size(); i++) {
        auto &f_i = mrcpp::get_func(xcInpVec, i);
        inpNodes += f_i.getNNodes();
        inpSize += f_i.getSizeNodes();
    }
    mrcpp::print::tree(3, "Preprocess input", inpNodes, inpSize, t_pre.elapsed());
    mrcpp::Timer t_eval;

    // divide nNodes into parts assigned to each MPI rank
    int nNodes = grid().size();
    int n_start = (mrcpp::mpi::wrk_rank * nNodes) / mrcpp::mpi::wrk_size;
    int n_end = ((mrcpp::mpi::wrk_rank + 1) * nNodes) / mrcpp::mpi::wrk_size;
    std::vector<Eigen::MatrixXd> ctrOutDataVec(n_end - n_start);
    mrcpp::FunctionTreeVector<3> ctrOutVec;
    if (mrcpp::mpi::wrk_size == 1) ctrOutVec = grid().generate(nOutCtr);

#pragma omp parallel
    {
#pragma omp for schedule(guided)
        for (int n = n_start; n < n_end; n++) {
            auto xcInpNodes = xc_utils::fetch_nodes(n, xcInpVec);// vector<mrcpp::FunctionNode<3> *>
            auto xcInpData = xc_utils::compress_nodes(xcInpNodes);// Eigen::MatrixXd
            auto xcOutData = functional().evaluate(xcInpData); // Eigen::MatrixXd
            auto ctrInpNodes = xc_utils::fetch_nodes(n, ctrInpVec);
            auto ctrInpData = xc_utils::compress_nodes(ctrInpNodes);// Eigen::MatrixXd
            auto ctrOutData = functional().contract(xcOutData, ctrInpData);// Eigen::MatrixXd .contract multiplies density and functional derivatives

            if (mrcpp::mpi::wrk_size > 1) {
                // store the results temporarily
                ctrOutDataVec[n - n_start] = std::move(ctrOutData);
            } else {
                // postprocess the results
                auto ctrOutNodes = xc_utils::fetch_nodes(n, ctrOutVec);
                xc_utils::expand_nodes(ctrOutNodes, ctrOutData);
            }
        }
    }

    // Input data is cleared before constructing the full output
    mrcpp::clear(xcInpVec, false);
    mrcpp::clear(ctrInpVec, false);

    if (mrcpp::mpi::wrk_size > 1) {
        // each MPI process has only a part of the results

        ctrOutVec = grid().generate(nOutCtr);
        mrcpp::BankAccount ctrOutBank; // to put the ctrOutDataVec;

        // note that mpi cannot run in multiple omp threads
        int size = nOutCtr * nCoefs;
        for (int n = n_start; n < n_end; n++) ctrOutBank.put_data(n, size, ctrOutDataVec[n - n_start].data());
        // fetch all nodes from bank and postprocess
        for (int n = 0; n < nNodes; n++) {
            Eigen::MatrixXd ctrOutData(nOutCtr, nCoefs);
            ctrOutBank.get_data(n, size, ctrOutData.data());
            auto ctrOutNodes = xc_utils::fetch_nodes(n, ctrOutVec);
            xc_utils::expand_nodes(ctrOutNodes, ctrOutData);
        }
    }

    // Reconstruct raw xcfun output functions
    /*
    for (auto i = 0; i < xcOutVec.size(); i++) {
        auto &f_i = mrcpp::get_func(xcOutVec, i);
        f_i.mwTransform(mrcpp::BottomUp);
        f_i.calcSquareNorm();
    }
    mrcpp::clear(xcOutVec, true);
    */

    // Reconstruct contracted output functions
    int ctrNodes = 0;
    int ctrSize = 0;
    for (auto i = 0; i < ctrOutVec.size(); i++) {
        auto &f_i = mrcpp::get_func(ctrOutVec, i);
        ctrNodes += f_i.getNNodes();
        ctrSize += f_i.getSizeNodes();
        f_i.mwTransform(mrcpp::BottomUp);
        f_i.calcSquareNorm();
    }
    mrcpp::print::tree(3, "Evaluate functional", ctrNodes, ctrSize, t_eval.elapsed());
    mrcpp::Timer t_post;
    auto potOutVec = functional().postprocess(ctrOutVec);
    mrcpp::clear(ctrOutVec, true);
    functional().clear();

    int outNodes = 0;
    int outSize = 0;

    for (auto i = 0; i < potOutVec.size(); i++) {
        auto &f_i = mrcpp::get_func(potOutVec, i);
        f_i.mwTransform(mrcpp::BottomUp);
        f_i.calcSquareNorm();
        // TODO? insert a crop
        outNodes += f_i.getNNodes();
        outSize += f_i.getSizeNodes();
    }

    mrcpp::print::tree(3, "Postprocess potential", outNodes, outSize, t_post.elapsed());
    return potOutVec;
    }
}

} // namespace mrdft
