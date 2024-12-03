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

    // Note: we do not need to make the complete energy density (PotVec[0]) for each mpi,
    //       instead we compute the energy directly.
    mrcpp::Timer t_post;
    int potvecSize = 2; // this size include PotVec[0], which is not used
    if (functional().isSpin()) potvecSize = 3;
    mrcpp::FunctionTreeVector<3> PotVec = grid().generate(potvecSize);
    int nNodes = grid().size();
    // parallelization of loop both with omp (pragma omp for) and
    // mpi (each mpi has a portion of the loop, defined by n_start and n_end)
    int n_start = (mrcpp::mpi::wrk_rank * nNodes) / mrcpp::mpi::wrk_size;
    int n_end = ((mrcpp::mpi::wrk_rank + 1) * nNodes) / mrcpp::mpi::wrk_size;
    DoubleVector XCenergy = DoubleVector::Zero(1);
    double sum = 0.0;
#pragma omp parallel
    {
#pragma omp for schedule(guided) reduction (+: sum)
        for (int n = n_start; n < n_end; n++) {
            vector<mrcpp::FunctionNode<3> *> xcNodes = xc_utils::fetch_nodes(n, PotVec);
            functional().makepot(inp, xcNodes);
            sum += xcNodes[0]->integrate();
        }
    }
    XCenergy[0] = sum;

    // each mpi only has part of the results. All send their results to bank and then fetch
    if(mrcpp::mpi::wrk_size > 1) {
        // sum up the energy contrbutions from all mpi
        mrcpp::mpi::allreduce_vector(XCenergy, mrcpp::mpi::comm_wrk);
        // send to bank, note that this cannot be done in the omp loop above,
        // because omp threads cannot use mpi
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
    this->functional().XCenergy = XCenergy[0];
    functional().clear();
    int outNodes = 0;
    int outSize = 0;
    for (int i = 1; i < PotVec.size(); i++) {
        mrcpp::FunctionTree<3> &f_i = mrcpp::get_func(PotVec, i);
        f_i.mwTransform(mrcpp::BottomUp);
        f_i.calcSquareNorm();
        outNodes += f_i.getNNodes();
        outSize += f_i.getSizeNodes();
    }
    mrcpp::print::tree(3, "Make potential", outNodes, outSize, t_post.elapsed());
    return PotVec;
}

} // namespace mrdft
