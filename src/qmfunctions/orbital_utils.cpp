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
#include <MRCPP/Timer>
#include <MRCPP/trees/FunctionNode.h>
#include <MRCPP/utils/details.h>

#include "parallel.h"
#include "utils/RRMaximizer.h"
#include "utils/math_utils.h"
#include "utils/print_utils.h"

#include "Orbital.h"
#include "OrbitalIterator.h"
#include "orbital_utils.h"
#include "qmfunction_utils.h"

using mrcpp::FunctionNode;
using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

namespace orbital {
ComplexMatrix localize(double prec, OrbitalVector &Phi, int spin);
ComplexMatrix calc_localization_matrix(double prec, OrbitalVector &Phi);
} // namespace orbital

/****************************************
 * Orbital related standalone functions *
 ****************************************/

/** @brief Compute <bra|ket> = int bra^\dag(r) * ket(r) dr.
 *
 *  Notice that the <bra| position is already complex conjugated.
 *  Alpha spin is orthogonal to beta spin, but paired orbitals are
 *  not necessarily orthogonal to alpha/beta orbitals.
 *
 */
ComplexDouble orbital::dot(Orbital bra, Orbital ket) {
    if ((bra.spin() == SPIN::Alpha) and (ket.spin() == SPIN::Beta)) return 0.0;
    if ((bra.spin() == SPIN::Beta) and (ket.spin() == SPIN::Alpha)) return 0.0;
    return qmfunction::dot(bra, ket);
}

/** @brief Compute the diagonal dot products <bra_i|ket_i>
 *
 * MPI: dot product is computed by the ket owner and the corresponding
 *      bra is communicated. The resulting vector is allreduced, and
 *      the foreign bra's are cleared.
 *
 */
ComplexVector orbital::dot(OrbitalVector &Bra, OrbitalVector &Ket) {
    if (Bra.size() != Ket.size()) MSG_ABORT("Size mismatch");

    int N = Bra.size();
    ComplexVector result = ComplexVector::Zero(N);
    for (int i = 0; i < N; i++) {
        // The bra is sent to the owner of the ket
        if (Bra[i].rankID() != Ket[i].rankID()) {
            int tag = 8765 + i;
            int src = Bra[i].rankID();
            int dst = Ket[i].rankID();
            if (mpi::my_orb(Bra[i])) mpi::send_function(Bra[i], dst, tag, mpi::comm_orb);
            if (mpi::my_orb(Ket[i])) mpi::recv_function(Bra[i], src, tag, mpi::comm_orb);
        }
        result[i] = orbital::dot(Bra[i], Ket[i]);
        if (not mpi::my_orb(Bra[i])) Bra[i].free(NUMBER::Total);
    }
    mpi::allreduce_vector(result, mpi::comm_orb);
    return result;
}

/** @brief Compute <bra|ket> = int |bra^\dag(r)| * |ket(r)| dr.
 *
 */
ComplexDouble orbital::node_norm_dot(Orbital bra, Orbital ket, bool exact) {
    if ((bra.spin() == SPIN::Alpha) and (ket.spin() == SPIN::Beta)) return 0.0;
    if ((bra.spin() == SPIN::Beta) and (ket.spin() == SPIN::Alpha)) return 0.0;
    return qmfunction::node_norm_dot(bra, ket, exact);
}

/** @brief Compare spin and occupation of two orbitals
 *
 *  Returns true if orbital parameters are the same.
 *
 */
bool orbital::compare(const Orbital &phi_a, const Orbital &phi_b) {
    bool comp = true;
    if (compare_occupation(phi_a, phi_b) < 0) {
        MSG_WARN("Different occupation");
        comp = false;
    }
    if (compare_spin(phi_a, phi_b) < 0) {
        MSG_WARN("Different spin");
        comp = false;
    }
    return comp;
}

/** @brief Compare occupation of two orbitals
 *
 *  Returns the common occupation if they match, -1 if they differ.
 *
 */
int orbital::compare_occupation(const Orbital &phi_a, const Orbital &phi_b) {
    int comp = -1;
    if (phi_a.occ() == phi_b.occ()) comp = phi_a.occ();
    return comp;
}

/** @brief Compare spin of two orbitals
 *
 *  Returns the common spin if they match, -1 if they differ.
 *
 */
int orbital::compare_spin(const Orbital &phi_a, const Orbital &phi_b) {
    int comp = -1;
    if (phi_a.spin() == phi_b.spin()) comp = phi_a.spin();
    return comp;
}

/** @brief Compare spin and occupation of two orbital vector
 *
 *  Returns true if orbital parameters are the same, orbital ordering
 *  NOT taken into account.
 *
 */
bool orbital::compare(const OrbitalVector &Phi_a, const OrbitalVector &Phi_b) {
    bool comp = true;
    if (orbital::size_alpha(Phi_a) != orbital::size_alpha(Phi_b)) {
        MSG_WARN("Different alpha occupancy");
        comp = false;
    }
    if (orbital::size_beta(Phi_a) != orbital::size_beta(Phi_b)) {
        MSG_WARN("Different beta occupancy");
        comp = false;
    }
    if (orbital::size_paired(Phi_a) != orbital::size_paired(Phi_b)) {
        MSG_WARN("Different paired occupancy");
        comp = false;
    }
    if (orbital::size_empty(Phi_a) != orbital::size_empty(Phi_b)) {
        MSG_WARN("Different empty occupancy");
        comp = false;
    }
    if (orbital::size_singly(Phi_a) != orbital::size_singly(Phi_b)) {
        MSG_WARN("Different single occupancy");
        comp = false;
    }
    if (orbital::size_doubly(Phi_a) != orbital::size_doubly(Phi_b)) {
        MSG_WARN("Different double occupancy");
        comp = false;
    }
    if (orbital::size_occupied(Phi_a) != orbital::size_occupied(Phi_b)) {
        MSG_WARN("Different total occupancy");
        comp = false;
    }

    for (auto &phi_a : Phi_a) {
        const mrcpp::MultiResolutionAnalysis<3> *mra_a{nullptr};
        if (phi_a.hasReal()) mra_a = &phi_a.real().getMRA();
        if (phi_a.hasImag()) mra_a = &phi_a.imag().getMRA();
        if (mra_a == nullptr) continue;
        for (auto &phi_b : Phi_b) {
            const mrcpp::MultiResolutionAnalysis<3> *mra_b{nullptr};
            if (phi_b.hasReal()) mra_b = &phi_a.real().getMRA();
            if (phi_b.hasImag()) mra_b = &phi_a.imag().getMRA();
            if (mra_b == nullptr) continue;
            if (*mra_a != *mra_b) {
                MSG_WARN("Different MRA");
                comp = false;
            }
        }
    }
    return comp;
}

/** @brief out_i = a*(inp_a)_i + b*(inp_b)_i
 *
 *  Component-wise addition of orbitals.
 *
 */
OrbitalVector orbital::add(ComplexDouble a, OrbitalVector &Phi_a, ComplexDouble b, OrbitalVector &Phi_b, double prec) {
    if (Phi_a.size() != Phi_b.size()) MSG_ERROR("Size mismatch");

    OrbitalVector out = orbital::param_copy(Phi_a);
    for (int i = 0; i < Phi_a.size(); i++) {
        if (Phi_a[i].rankID() != Phi_b[i].rankID()) MSG_ABORT("MPI rank mismatch");
        qmfunction::add(out[i], a, Phi_a[i], b, Phi_b[i], prec);
    }
    return out;
}

/** @brief Orbital transformation out_j = sum_i inp_i*U_ij
 *
 * NOTE: OrbitalVector is considered a ROW vector, so rotation
 *       means matrix multiplication from the right
 *
 * MPI: Rank distribution of output vector is the same as input vector
 *
 */
OrbitalVector orbital::rotate(OrbitalVector &Phi, const ComplexMatrix &U, double prec) {

    // The principle of this routine is that nodes are rotated one by one using matrix multiplication.
    // The routine does avoid when possible to move data, but uses pointers and indices manipulation.
    // MPI version does not use OMP yet, Serial version uses OMP

    auto out = orbital::param_copy(Phi);
    int N = Phi.size();

    // 1) make union tree without coefficients
    mrcpp::FunctionTree<3> refTree(*MRA);
    mpi::allreduce_Tree_noCoeff(refTree, Phi, mpi::comm_orb);

    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();
    std::vector<double> scalefac_ref;
    std::vector<double *> coeffVec_ref; // not used!
    std::vector<int> indexVec_ref;      // serialIx of the nodes
    std::vector<int> parindexVec_ref;   // serialIx of the parent nodes
    int max_ix;
    // get a list of all nodes in union tree, identified by their serialIx indices
    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac_ref, max_ix, refTree);
    int max_n = indexVec_ref.size();

    // 2) We work with real numbers only. Make real blocks for U matrix
    bool UhasReal = false;
    bool UhasImag = false;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (std::abs(U(i, j).real()) > mrcpp::MachineZero) UhasReal = true;
            if (std::abs(U(i, j).imag()) > mrcpp::MachineZero) UhasImag = true;
        }
    }

    IntVector PsihasReIm = IntVector::Zero(2);
    for (int i = 0; i < N; i++) {
        if (!mpi::my_orb(Phi[i])) continue;
        PsihasReIm[0] = (Phi[i].hasReal()) ? 1 : 0;
        PsihasReIm[1] = (Phi[i].hasImag()) ? 1 : 0;
    }
    mpi::allreduce_vector(PsihasReIm, mpi::comm_orb);
    if (not PsihasReIm[0] and not PsihasReIm[1]) {
        OrbitalVector out = orbital::param_copy(Phi);
        return out;
    }

    bool makeReal = (UhasReal and PsihasReIm[0]) or (UhasImag and PsihasReIm[1]);
    bool makeImag = (UhasReal and PsihasReIm[1]) or (UhasImag and PsihasReIm[0]);

    if (not makeReal and not makeImag) {
        OrbitalVector out = orbital::param_copy(Phi);
        return out;
    }

    int Neff = N;               // effective number of orbitals
    if (makeImag) Neff = 2 * N; // Imag and Real treated independently. We always use real part of U

    IntVector conjMat = IntVector::Zero(Neff);
    for (int i = 0; i < Neff; i++) {
        if (!mpi::my_orb(Phi[i % N])) continue;
        conjMat[i] = (Phi[i % N].conjugate()) ? -1 : 1;
    }
    mpi::allreduce_vector(conjMat, mpi::comm_orb);

    // we make a real matrix = U,  but organized as one or four real blocks
    // out_r = U_rr*in_r - U_ir*in_i*conjMat
    // out_i = U_ri*in_r - U_ii*in_i*conjMat
    // the first index of U is the one used on input Phi
    DoubleMatrix Ureal(Neff, Neff); // four blocks, for rr ri ir ii
    for (int i = 0; i < Neff; i++) {
        for (int j = 0; j < Neff; j++) {
            double sign = 1.0;
            if (i < N and j < N) {
                // real U applied on real Phi
                Ureal(i, j) = U.real()(i % N, j % N);
            } else if (i >= N and j >= N) {
                // real U applied on imag Phi
                Ureal(i, j) = conjMat[i] * U.real()(i % N, j % N);
            } else if (i < N and j >= N) {
                // imag U applied on real Phi
                Ureal(i, j) = U.imag()(i % N, j % N);
            } else {
                // imag U applied on imag Phi
                Ureal(i, j) = -1.0 * conjMat[i] * U.imag()(i % N, j % N);
            }
        }
    }

    // 3) In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank

    bool serial = mpi::orb_size == 1; // flag for serial/MPI switch
    BankAccount nodesPhi;             // to put the original nodes
    BankAccount nodesRotated;         // to put the rotated nodes

    // used for serial only:
    std::vector<std::vector<double *>> coeffVec(Neff);
    std::vector<std::vector<int>> indexVec(Neff);   // serialIx of the nodes
    std::map<int, std::vector<int>> node2orbVec;    // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2node(Neff); // for a given orbital and a given node, gives the node index in the
                                                    // orbital given the node index in the reference tree

    if (serial) {

        // make list of all coefficients (coeffVec), and their reference indices (indexVec)
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<double> scalefac;
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            if (Phi[j].hasReal()) {
                Phi[j].real().makeCoeffVector(coeffVec[j], indexVec[j], parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec[j]) {
                    orb2node[j][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j);
                }
            }
            if (Phi[j].hasImag()) {
                Phi[j].imag().makeCoeffVector(coeffVec[j + N], indexVec[j + N], parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec[j + N]) {
                    orb2node[j + N][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j + N);
                }
            }
        }

    } else { // MPI case

        // send own nodes to bank, identifying them through the serialIx of refTree
        save_nodes(Phi, refTree, nodesPhi);
        mrchem::mpi::barrier(mrchem::mpi::comm_orb); // required for now, as the blockdata functionality has no queue yet.
    }

    // 4) rotate all the nodes
    IntMatrix split_serial;                             // in the serial case all split are store in one array
    std::vector<double> split(Neff, -1.0);              // which orbitals need splitting (at a given node). For now double for compatibilty with bank
    std::vector<double> needsplit(Neff, 1.0);           // which orbitals need splitting
    std::vector<std::vector<double *>> coeffpVec(Neff); // to put pointers to the rotated coefficient for each orbital
    std::vector<std::map<int, int>> ix2coef(Neff);      // to find the index in for example rotCoeffVec[] corresponding to a serialIx
    int csize;                                          // size of the current coefficients (different for roots and branches)
    std::vector<DoubleMatrix> rotatedCoeffVec;          // just to ensure that the data from rotatedCoeff is not deleted, since we point to it.
    // j indices are for unrotated orbitals, i indices are for rotated orbitals
    if (serial) {
        std::map<int, int> ix2coef_ref;   // to find the index n corresponding to a serialIx
        split_serial.resize(Neff, max_n); // not use in the MPI case
        for (int n = 0; n < max_n; n++) {
            int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree
            ix2coef_ref[node_ix] = n;
            for (int i = 0; i < Neff; i++) split_serial(i, n) = 1;
        }

        std::vector<int> nodeReady(max_n, 0); // To indicate to OMP threads that the parent is ready (for splits)

        // assumes the nodes are ordered such that parent are treated before children. BFS or DFS ok.
        // NB: the n must be traversed approximately in right order: Thread n may have to wait until som other preceding
        // n is finished.
#pragma omp parallel for schedule(dynamic)
        for (int n = 0; n < max_n; n++) {
            int csize;
            int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree
            // 4a) make a dense contiguous matrix with the coefficient from all the orbitals using node n
            std::vector<int> orbjVec; // to remember which orbital correspond to each orbVec.size();
            if (node2orbVec[node_ix].size() <= 0) continue;
            csize = sizecoeffW;
            if (parindexVec_ref[n] < 0) csize = sizecoeff; // for root nodes we include scaling coeff

            int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
            if (parindexVec_ref[n] < 0) shift = 0;
            DoubleMatrix coeffBlock(csize, node2orbVec[node_ix].size());
            for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                int orb_node_ix = orb2node[j][node_ix];
                for (int k = 0; k < csize; k++) coeffBlock(k, orbjVec.size()) = coeffVec[j][orb_node_ix][k + shift];
                orbjVec.push_back(j);
            }

            // 4b) make a list of rotated orbitals needed for this node
            // OMP must wait until parent is ready
            while (parindexVec_ref[n] >= 0 and nodeReady[ix2coef_ref[parindexVec_ref[n]]] == 0) {
#pragma omp flush
            };

            std::vector<int> orbiVec;
            for (int i = 0; i < Neff; i++) { // loop over all rotated orbitals
                if (not makeReal and i < N) continue;
                if (not makeImag and i >= N) continue;
                if (parindexVec_ref[n] >= 0 and split_serial(i, ix2coef_ref[parindexVec_ref[n]]) == 0) continue; // parent node has too small wavelets
                orbiVec.push_back(i);
            }

            // 4c) rotate this node
            DoubleMatrix Un(orbjVec.size(), orbiVec.size()); // chunk of U, with reorganized indices
            for (int i = 0; i < orbiVec.size(); i++) {       // loop over rotated orbitals
                for (int j = 0; j < orbjVec.size(); j++) { Un(j, i) = Ureal(orbjVec[j], orbiVec[i]); }
            }
            DoubleMatrix rotatedCoeff(csize, orbiVec.size());
            // HERE IT HAPPENS!
            rotatedCoeff.noalias() = coeffBlock * Un; // Matrix mutiplication

            // 4d) store and make rotated node pointers
            // for now we allocate in buffer, in future could be directly allocated in the final trees
            double thres = prec * prec * scalefac_ref[n] * scalefac_ref[n];
            // make all norms:
            for (int i = 0; i < orbiVec.size(); i++) {
                // check if parent must be split
                if (parindexVec_ref[n] == -1 or split_serial(orbiVec[i], ix2coef_ref[parindexVec_ref[n]])) {
                    // mark this node for this orbital for later split
#pragma omp critical
                    {
                        ix2coef[orbiVec[i]][node_ix] = coeffpVec[orbiVec[i]].size();
                        coeffpVec[orbiVec[i]].push_back(&(rotatedCoeff(0, i))); // list of coefficient pointers
                    }
                    // check norms for split
                    double wnorm = 0.0; // rotatedCoeff(k, i) is already in cache here
                    int kstart = 0;
                    if (parindexVec_ref[n] < 0) kstart = sizecoeff - sizecoeffW; // do not include scaling, even for roots
                    for (int k = kstart; k < csize; k++) wnorm += rotatedCoeff(k, i) * rotatedCoeff(k, i);
                    if (thres < wnorm or prec < 0)
                        split_serial(orbiVec[i], n) = 1;
                    else
                        split_serial(orbiVec[i], n) = 0;
                } else {
                    ix2coef[orbiVec[i]][node_ix] = max_n + 1; // should not be used
                    split_serial(orbiVec[i], n) = 0;          // do not split if parent does not need to be split
                }
            }
            nodeReady[n] = 1;
#pragma omp critical
            {
                // this ensures that rotatedCoeff is not deleted, when getting out of scope
                rotatedCoeffVec.push_back(std::move(rotatedCoeff));
            }
        }

    } else { // MPI case

        // TODO? rotate in bank, so that we do not get and put. Requires clever handling of splits.
        std::vector<double> split(Neff, -1.0);    // which orbitals need splitting (at a given node). For now double for compatibilty with bank
        std::vector<double> needsplit(Neff, 1.0); // which orbitals need splitting
        BankAccount nodeSplits;
        mrchem::mpi::barrier(mrchem::mpi::comm_orb); // required for now, as the blockdata functionality has no queue yet.

        DoubleMatrix coeffBlock(sizecoeff, Neff);
        max_ix++; // largest node index + 1. to store rotated orbitals with different id
        TaskManager tasks(max_n);
        for (int nn = 0; nn < max_n; nn++) {
            int n = tasks.next_task();
            if (n < 0) break;
            // if (n % mpi::orb_size != mpi::orb_rank) continue; // could use any partitioning
            double thres = prec * prec * scalefac_ref[n] * scalefac_ref[n];
            // 4a) make list of orbitals that should split the parent node, i.e. include this node
            int parentid = parindexVec_ref[n];
            if (parentid == -1) {
                // root node, split if output needed
                for (int i = 0; i < N; i++) {
                    if (makeReal)
                        split[i] = 1.0;
                    else
                        split[i] = -1.0;
                }
                for (int i = N; i < Neff; i++) {
                    if (makeImag)
                        split[i] = 1.0;
                    else
                        split[i] = -1.0;
                }
                csize = sizecoeff;
            } else {
                // note that it will wait until data is available
                nodeSplits.get_data(parentid, Neff, split.data());
                csize = sizecoeffW;
            }
            std::vector<int> orbiVec;
            std::vector<int> orbjVec;
            for (int i = 0; i < Neff; i++) {  // loop over rotated orbitals
                if (split[i] < 0.0) continue; // parent node has too small wavelets
                orbiVec.push_back(i);
            }

            // 4b) rotate this node
            DoubleMatrix coeffBlock(csize, Neff); // largest possible used size
            nodesPhi.get_nodeblock(indexVec_ref[n], coeffBlock.data(), orbjVec);
            coeffBlock.conservativeResize(Eigen::NoChange, orbjVec.size()); // keep only used part

            // chunk of U, with reorganized indices and separate blocks for real and imag:
            DoubleMatrix Un(orbjVec.size(), orbiVec.size());
            DoubleMatrix rotatedCoeff(csize, orbiVec.size());

            for (int i = 0; i < orbiVec.size(); i++) {     // loop over included rotated real and imag part of orbitals
                for (int j = 0; j < orbjVec.size(); j++) { // loop over input orbital, possibly imaginary parts
                    Un(j, i) = Ureal(orbjVec[j], orbiVec[i]);
                }
            }

            // HERE IT HAPPENS
            rotatedCoeff.noalias() = coeffBlock * Un; // Matrix mutiplication

            // 3c) find which orbitals need to further refine this node, and store rotated node (after each other while
            // in cache).
            for (int i = 0; i < orbiVec.size(); i++) { // loop over rotated orbitals
                needsplit[orbiVec[i]] = -1.0;          // default, do not split
                // check if this node/orbital needs further refinement
                double wnorm = 0.0;
                int kwstart = csize - sizecoeffW; // do not include scaling
                for (int k = kwstart; k < csize; k++) wnorm += rotatedCoeff.col(i)[k] * rotatedCoeff.col(i)[k];
                if (thres < wnorm or prec < 0) needsplit[orbiVec[i]] = 1.0;
                nodesRotated.put_nodedata(orbiVec[i], indexVec_ref[n] + max_ix, csize, rotatedCoeff.col(i).data());
            }
            nodeSplits.put_data(indexVec_ref[n], Neff, needsplit.data());
        }
        mpi::barrier(mpi::comm_orb); // wait until all rotated nodes are ready
    }

    // 5) reconstruct trees using rotated nodes.

    // only serial case can use OMP, because MPI cannot be used by threads
    if (serial) {
        // OMP parallelized, but does not scale well, because the total memory bandwidth is a bottleneck. (the main
        // operation is writing the coefficient into the tree)

#pragma omp parallel for schedule(static)
        for (int j = 0; j < Neff; j++) {
            if (j < N) {
                out[j].alloc(NUMBER::Real);
                out[j].real().makeTreefromCoeff(refTree, coeffpVec[j], ix2coef[j], prec);
            } else {
                out[j % N].alloc(NUMBER::Imag);
                out[j % N].imag().makeTreefromCoeff(refTree, coeffpVec[j], ix2coef[j], prec);
            }
        }

    } else { // MPI case

        for (int j = 0; j < Neff; j++) {
            if (not mpi::my_orb(out[j % N])) continue;
            // traverse possible nodes, and stop descending when norm is zero (leaf in out[j])
            std::vector<double *> coeffpVec; //
            std::map<int, int> ix2coef;      // to find the index in coeffVec[] corresponding to a serialIx
            int ix = 0;
            std::vector<double *> pointerstodelete; // list of temporary arrays to clean up
            for (int ibank = 0; ibank < mpi::bank_size; ibank++) {
                std::vector<int> nodeidVec;
                double *dataVec; // will be allocated by bank
                nodesRotated.get_orbblock(j, dataVec, nodeidVec, ibank);
                if (nodeidVec.size() > 0) pointerstodelete.push_back(dataVec);
                int shift = 0;
                for (int n = 0; n < nodeidVec.size(); n++) {
                    if (nodeidVec[n] - max_ix < 0) { std::cout << n << " ERROR " << nodeidVec[n] << " " << max_ix << std::endl; }
                    assert(nodeidVec[n] - max_ix >= 0); // unrotated nodes have been deleted
                    if (ix2coef.count(nodeidVec[n] - max_ix) != 0) { std::cout << n << " ERRORix2 " << ix2coef.count(nodeidVec[n]) << " " << max_ix << std::endl; }
                    assert(ix2coef.count(nodeidVec[n] - max_ix) == 0); // each nodeid treated once
                    ix2coef[nodeidVec[n] - max_ix] = ix++;
                    csize = sizecoeffW;
                    if (parindexVec_ref[nodeidVec[n] - max_ix] < 0) csize = sizecoeff;
                    coeffpVec.push_back(&dataVec[shift]); // list of coeff pointers
                    shift += csize;
                }
            }
            if (j < N) {
                // Real part
                out[j].alloc(NUMBER::Real);
                out[j].real().makeTreefromCoeff(refTree, coeffpVec, ix2coef, prec);
            } else {
                // Imag part
                out[j % N].alloc(NUMBER::Imag);
                out[j % N].imag().makeTreefromCoeff(refTree, coeffpVec, ix2coef, prec);
            }
            for (double *p : pointerstodelete) delete[] p;
            pointerstodelete.clear();
        }
    }
    return out;
}

/** @brief Save all nodes in bank; identify them using serialIx from refTree
 * shift is a shift applied in the id
 */
void orbital::save_nodes(OrbitalVector Phi, mrcpp::FunctionTree<3> &refTree, BankAccount &account) {
    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();
    int max_nNodes = refTree.getNNodes();
    std::vector<double *> coeffVec;
    std::vector<double> scalefac;
    std::vector<int> indexVec;    // SerialIx of the node in refOrb
    std::vector<int> parindexVec; // SerialIx of the parent node
    int N = Phi.size();
    int max_ix;
    for (int j = 0; j < N; j++) {
        if (not mpi::my_orb(Phi[j])) continue;
        // make vector with all coef address and their index in the union grid
        if (Phi[j].hasReal()) {
            Phi[j].real().makeCoeffVector(coeffVec, indexVec, parindexVec, scalefac, max_ix, refTree);
            int max_n = indexVec.size();
            // send node coefs from Phi[j] to bank
            // except for the root nodes, only wavelets are sent
            for (int i = 0; i < max_n; i++) {
                if (indexVec[i] < 0) continue; // nodes that are not in refOrb
                int csize = sizecoeffW;
                if (parindexVec[i] < 0) csize = sizecoeff;
                account.put_nodedata(j, indexVec[i], csize, &(coeffVec[i][sizecoeff - csize]));
            }
        }
        // Imaginary parts are considered as orbitals with an orbid shifted by N
        if (Phi[j].hasImag()) {
            Phi[j].imag().makeCoeffVector(coeffVec, indexVec, parindexVec, scalefac, max_ix, refTree);
            int max_n = indexVec.size();
            // send node coefs from Phi[j] to bank
            for (int i = 0; i < max_n; i++) {
                if (indexVec[i] < 0) continue; // nodes that are not in refOrb
                // NB: the identifier (indexVec[i]) must be shifted for not colliding with the nodes from the real part
                int csize = sizecoeffW;
                if (parindexVec[i] < 0) csize = sizecoeff;
                account.put_nodedata(j + N, indexVec[i], csize, &(coeffVec[i][sizecoeff - csize]));
            }
        }
    }
}

/** @brief Deep copy
 *
 * New orbitals are constructed as deep copies of the input set.
 *
 */
OrbitalVector orbital::deep_copy(OrbitalVector &Phi) {
    OrbitalVector out;
    for (auto &i : Phi) {
        Orbital out_i = i.paramCopy();
        if (mpi::my_orb(out_i)) qmfunction::deep_copy(out_i, i);
        out.push_back(out_i);
    }
    return out;
}

/** @brief Parameter copy
 *
 * New orbitals are constructed as parameter copies of the input set.
 *
 */
OrbitalVector orbital::param_copy(const OrbitalVector &Phi) {
    OrbitalVector out;
    for (const auto &i : Phi) {
        Orbital out_i = i.paramCopy();
        out.push_back(out_i);
    }
    return out;
}

/** @brief Adjoin two vectors
 *
 * The orbitals of the input vector are appended to
 * (*this) vector, the ownership is transferred. Leaves
 * the input vector empty.
 *
 */
OrbitalVector orbital::adjoin(OrbitalVector &Phi_a, OrbitalVector &Phi_b) {
    OrbitalVector out;
    for (auto &phi : Phi_a) out.push_back(phi);
    for (auto &phi : Phi_b) out.push_back(phi);
    Phi_a.clear();
    Phi_b.clear();
    return out;
}

/** @brief Disjoin vector in two parts
 *
 * All orbitals of a particular spin is collected in a new vector
 * and returned. These orbitals are removed from (*this) vector,
 * and the ownership is transferred.
 *
 */
OrbitalVector orbital::disjoin(OrbitalVector &Phi, int spin) {
    OrbitalVector out;
    OrbitalVector tmp;
    for (auto &i : Phi) {
        if (i.spin() == spin) {
            out.push_back(i);
        } else {
            tmp.push_back(i);
        }
    }
    Phi.clear();
    Phi = tmp;
    return out;
}

/** @brief Write orbitals to disk
 *
 * @param Phi: orbitals to save
 * @param file: file name prefix
 * @param spin: type of orbitals to save, negative means all orbitals
 *
 * The given file name (e.g. "phi") will be appended with orbital number ("phi_0").
 * Produces separate files for meta data ("phi_0.meta"), real ("phi_0_re.tree") and
 * imaginary ("phi_0_im.tree") parts. If a particular spin is given, the file name
 * will get an extra "_p", "_a" or "_b" suffix. Negative spin means that all
 * orbitals in the vector are saved, and no suffix is added.
 */
void orbital::save_orbitals(OrbitalVector &Phi, const std::string &file, int spin) {
    Timer t_tot;
    std::string spin_str = "All";
    if (spin == SPIN::Paired) spin_str = "Paired";
    if (spin == SPIN::Alpha) spin_str = "Alpha";
    if (spin == SPIN::Beta) spin_str = "Beta";
    mrcpp::print::header(2, "Writing orbitals");
    print_utils::text(2, "File name", file);
    print_utils::text(2, "Spin", spin_str);
    mrcpp::print::separator(2, '-');

    auto n = 0;
    for (int i = 0; i < Phi.size(); i++) {
        if ((Phi[i].spin() == spin) or (spin < 0)) {
            Timer t1;
            std::stringstream orbname;
            orbname << file << "_idx_" << n;
            if (mpi::my_orb(Phi[i])) Phi[i].saveOrbital(orbname.str());
            print_utils::qmfunction(2, "'" + orbname.str() + "'", Phi[i], t1);
            n++;
        }
    }
    mrcpp::print::footer(2, t_tot, 2);
}

/** @brief Read orbitals from disk
 *
 * @param file: file name prefix
 * @param n_orbs: number of orbitals to read
 *
 * The given file name (e.g. "phi") will be appended with orbital number ("phi_0").
 * Reads separate files for meta data ("phi_0.meta"), real ("phi_0_re.tree") and
 * imaginary ("phi_0_im.tree") parts. Negative n_orbs means that all orbitals matching
 * the prefix name will be read.
 */
OrbitalVector orbital::load_orbitals(const std::string &file, int n_orbs) {
    Timer t_tot;
    mrcpp::print::header(2, "Reading orbitals");
    print_utils::text(2, "File name", file);
    mrcpp::print::separator(2, '-');
    OrbitalVector Phi;
    for (int i = 0; true; i++) {
        if (n_orbs > 0 and i >= n_orbs) break;
        Timer t1;
        Orbital phi_i;
        std::stringstream orbname;
        orbname << file << "_idx_" << i;
        phi_i.loadOrbital(orbname.str());
        phi_i.setRankID(mpi::orb_rank);
        if (phi_i.hasReal() or phi_i.hasImag()) {
            phi_i.setRankID(i % mpi::orb_size);
            Phi.push_back(phi_i);
            print_utils::qmfunction(2, "'" + orbname.str() + "'", phi_i, t1);
            if (not mpi::my_orb(phi_i)) phi_i.free(NUMBER::Total);
        } else {
            break;
        }
    }
    mrcpp::print::footer(2, t_tot, 2);
    return Phi;
}

/** @brief Normalize single orbital. Private function. */
void orbital::normalize(Orbital &phi) {
    phi.rescale(1.0 / phi.norm());
}

/** @brief Normalize all orbitals in the set */
void orbital::normalize(OrbitalVector &Phi) {
    mpi::free_foreign(Phi);
    for (auto &phi_i : Phi)
        if (mpi::my_orb(phi_i)) orbital::normalize(phi_i);
}

/** @brief In place orthogonalize against inp. Private function. */
void orbital::orthogonalize(double prec, Orbital &phi, Orbital psi) {
    ComplexDouble overlap = orbital::dot(psi, phi);
    double sq_norm = psi.squaredNorm();
    if (std::abs(overlap) > prec) phi.add(-1.0 * overlap / sq_norm, psi);
}

/** @brief Gram-Schmidt orthogonalize orbitals within the set */
void orbital::orthogonalize(double prec, OrbitalVector &Phi) {
    mpi::free_foreign(Phi);
    for (int i = 0; i < Phi.size(); i++) {
        for (int j = 0; j < i; j++) {
            int tag = 7632 * i + j;
            int src = Phi[j].rankID();
            int dst = Phi[i].rankID();
            if (src == dst) {
                if (mpi::my_orb(Phi[i])) orbital::orthogonalize(prec / Phi.size(), Phi[i], Phi[j]);
            } else {
                if (mpi::my_orb(Phi[i])) {
                    mpi::recv_function(Phi[j], src, tag, mpi::comm_orb);
                    orbital::orthogonalize(prec / Phi.size(), Phi[i], Phi[j]);
                    Phi[j].free(NUMBER::Total);
                }
                if (mpi::my_orb(Phi[j])) mpi::send_function(Phi[j], dst, tag, mpi::comm_orb);
            }
        }
    }
}

/** @brief Orthogonalize the Phi orbital against all orbitals in Psi */
void orbital::orthogonalize(double prec, OrbitalVector &Phi, OrbitalVector &Psi) {
    // Get all output orbitals belonging to this MPI
    OrbitalChunk myPhi = mpi::get_my_chunk(Phi);

    // Orthogonalize MY orbitals with ALL input orbitals
    OrbitalIterator iter(Psi, false);
    while (iter.next()) {
        for (int i = 0; i < iter.get_size(); i++) {
            Orbital &psi_i = iter.orbital(i);
            for (auto &j : myPhi) {
                Orbital &phi_j = std::get<1>(j);
                orbital::orthogonalize(prec / Psi.size(), phi_j, psi_i);
            }
        }
    }
}

/** @brief Orbital transformation out_j = sum_i inp_i*U_ij
 *
 * NOTE: OrbitalVector is considered a ROW vector, so rotation
 *       means matrix multiplication from the right
 *
 * MPI: Rank distribution of output vector is the same as input vector
 *
 */
ComplexMatrix orbital::calc_overlap_matrix(OrbitalVector &BraKet) {

    // TODO: spin separate in block?
    int N = BraKet.size();
    ComplexMatrix S = ComplexMatrix::Zero(N, N);
    DoubleMatrix Sreal = DoubleMatrix::Zero(2 * N, 2 * N); // same as S, but stored as 4 blocks, rr,ri,ir,ii

    // 1) make union tree without coefficients
    mrcpp::FunctionTree<3> refTree(*MRA);
    mpi::allreduce_Tree_noCoeff(refTree, BraKet, mpi::comm_orb);

    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();

    // get a list of all nodes in union grid, as defined by their indices
    std::vector<double> scalefac;
    std::vector<double *> coeffVec_ref;
    std::vector<int> indexVec_ref;    // serialIx of the nodes
    std::vector<int> parindexVec_ref; // serialIx of the parent nodes
    int max_ix;                       // largest index value (not used here)

    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac, max_ix, refTree);
    int max_n = indexVec_ref.size();

    // only used for serial case:
    std::vector<std::vector<double *>> coeffVec(2 * N);
    std::map<int, std::vector<int>> node2orbVec;     // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2node(2 * N); // for a given orbital and a given node, gives the node index in
                                                     // the orbital given the node index in the reference tree

    bool serial = mpi::orb_size == 1; // flag for serial/MPI switch
    BankAccount nodesBraKet;

    // In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank
    if (serial) {
        // 2) make list of all coefficients, and their reference indices
        // for different orbitals, indexVec will give the same index for the same node in space
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<int> indexVec;    // serialIx of the nodes
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            if (BraKet[j].hasReal()) {
                BraKet[j].real().makeCoeffVector(coeffVec[j], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2node[j][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j);
                }
            }
            if (BraKet[j].hasImag()) {
                BraKet[j].imag().makeCoeffVector(coeffVec[j + N], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2node[j + N][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j + N);
                }
            }
        }
    } else { // MPI case
        // 2) send own nodes to bank, identifying them through the serialIx of refTree
        save_nodes(BraKet, refTree, nodesBraKet);
        mrchem::mpi::barrier(mrchem::mpi::comm_orb); // wait until everything is stored before fetching!
    }

    // 3) make dot product for all the nodes and accumulate into S

    int ibank = 0;
#pragma omp parallel for schedule(dynamic) if (serial)
    for (int n = 0; n < max_n; n++) {
        if (n % mpi::orb_size != mpi::orb_rank) continue;
        int csize;
        int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree
        std::vector<int> orbVec;       // identifies which orbitals use this node
        if (serial and node2orbVec[node_ix].size() <= 0) continue;
        if (parindexVec_ref[n] < 0)
            csize = sizecoeff;
        else
            csize = sizecoeffW;

        // In the serial case we copy the coeff coeffBlock. In the mpi case coeffBlock is provided by the bank
        if (serial) {
            int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
            if (parindexVec_ref[n] < 0) shift = 0;
            DoubleMatrix coeffBlock(csize, node2orbVec[node_ix].size());
            for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                int orb_node_ix = orb2node[j][node_ix];
                for (int k = 0; k < csize; k++) coeffBlock(k, orbVec.size()) = coeffVec[j][orb_node_ix][k + shift];
                orbVec.push_back(j);
            }
            if (orbVec.size() > 0) {
                DoubleMatrix S_temp(orbVec.size(), orbVec.size());
                S_temp.noalias() = coeffBlock.transpose() * coeffBlock;
                for (int i = 0; i < orbVec.size(); i++) {
                    for (int j = 0; j < orbVec.size(); j++) {
                        if (BraKet[orbVec[i] % N].spin() == SPIN::Alpha and BraKet[orbVec[j] % N].spin() == SPIN::Beta) continue;
                        if (BraKet[orbVec[i] % N].spin() == SPIN::Beta and BraKet[orbVec[j] % N].spin() == SPIN::Alpha) continue;
                        double &Srealij = Sreal(orbVec[i], orbVec[j]);
                        double &Stempij = S_temp(i, j);
#pragma omp atomic
                        Srealij += Stempij;
                    }
                }
            }
        } else { // MPI case
            DoubleMatrix coeffBlock(csize, 2 * N);
            nodesBraKet.get_nodeblock(indexVec_ref[n], coeffBlock.data(), orbVec);

            if (orbVec.size() > 0) {
                DoubleMatrix S_temp(orbVec.size(), orbVec.size());
                coeffBlock.conservativeResize(Eigen::NoChange, orbVec.size());
                S_temp.noalias() = coeffBlock.transpose() * coeffBlock;
                for (int i = 0; i < orbVec.size(); i++) {
                    for (int j = 0; j < orbVec.size(); j++) {
                        if (BraKet[orbVec[i] % N].spin() == SPIN::Alpha and BraKet[orbVec[j] % N].spin() == SPIN::Beta) continue;
                        if (BraKet[orbVec[i] % N].spin() == SPIN::Beta and BraKet[orbVec[j] % N].spin() == SPIN::Alpha) continue;
                        Sreal(orbVec[i], orbVec[j]) += S_temp(i, j);
                    }
                }
            }
        }
    }

    IntVector conjMat = IntVector::Zero(N);
    for (int i = 0; i < N; i++) {
        if (!mpi::my_orb(BraKet[i])) continue;
        conjMat[i] = (BraKet[i].conjugate()) ? -1 : 1;
    }
    mpi::allreduce_vector(conjMat, mpi::comm_orb);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            S.real()(i, j) = Sreal(i, j) + conjMat[i] * conjMat[j] * Sreal(i + N, j + N);
            S.imag()(i, j) = conjMat[j] * Sreal(i, j + N) - conjMat[i] * Sreal(i + N, j);
            if (i != j) S(j, i) = std::conj(S(i, j)); // ensure exact symmetri
        }
    }

    // Assumes linearity: result is sum of all nodes contributions
    mpi::allreduce_matrix(S, mpi::comm_orb);

    return S;
}

/** @brief Compute the overlap matrix S_ij = <bra_i|ket_j>
 *
 */
ComplexMatrix orbital::calc_overlap_matrix(OrbitalVector &Bra, OrbitalVector &Ket) {

    // TODO: spin separate in block?
    int N = Bra.size();
    int M = Ket.size();
    ComplexMatrix S = ComplexMatrix::Zero(N, M);
    DoubleMatrix Sreal = DoubleMatrix::Zero(2 * N, 2 * M); // same as S, but stored as 4 blocks, rr,ri,ir,ii

    // 1) make union tree without coefficients for Bra (supposed smallest)
    mrcpp::FunctionTree<3> refTree(*MRA);
    mpi::allreduce_Tree_noCoeff(refTree, Bra, mpi::comm_orb);
    // note that Ket is not part of union grid: if a node is in ket but not in Bra, the dot product is zero.

    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();

    // get a list of all nodes in union grid, as defined by their indices
    std::vector<double *> coeffVec_ref;
    std::vector<int> indexVec_ref;    // serialIx of the nodes
    std::vector<int> parindexVec_ref; // serialIx of the parent nodes
    std::vector<double> scalefac;
    int max_ix;

    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac, max_ix, refTree);
    int max_n = indexVec_ref.size();
    max_ix++;

    bool serial = mpi::orb_size == 1; // flag for serial/MPI switch

    // only used for serial case:
    std::vector<std::vector<double *>> coeffVecBra(2 * N);
    std::map<int, std::vector<int>> node2orbVecBra;     // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2nodeBra(2 * N); // for a given orbital and a given node, gives the node index in
                                                        // the orbital given the node index in the reference tree
    std::vector<std::vector<double *>> coeffVecKet(2 * M);
    std::map<int, std::vector<int>> node2orbVecKet;     // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2nodeKet(2 * M); // for a given orbital and a given node, gives the node index in
                                                        // the orbital given the node index in the reference tree
    BankAccount nodesBra;
    BankAccount nodesKet;

    // In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank
    if (serial) {
        // 2) make list of all coefficients, and their reference indices
        // for different orbitals, indexVec will give the same index for the same node in space
        // TODO? : do not copy coefficients, but use directly the pointers
        // could OMP parallelize, but is fast anyway
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<int> indexVec;    // serialIx of the nodes
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            if (Bra[j].hasReal()) {
                Bra[j].real().makeCoeffVector(coeffVecBra[j], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2nodeBra[j][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVecBra[ix].push_back(j);
                }
            }
            if (Bra[j].hasImag()) {
                Bra[j].imag().makeCoeffVector(coeffVecBra[j + N], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2nodeBra[j + N][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVecBra[ix].push_back(j + N);
                }
            }
        }
        for (int j = 0; j < M; j++) {
            if (Ket[j].hasReal()) {
                Ket[j].real().makeCoeffVector(coeffVecKet[j], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2nodeKet[j][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVecKet[ix].push_back(j);
                }
            }
            if (Ket[j].hasImag()) {
                Ket[j].imag().makeCoeffVector(coeffVecKet[j + M], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2nodeKet[j + M][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVecKet[ix].push_back(j + M);
                }
            }
        }

    } else { // MPI case

        // 2) send own nodes to bank, identifying them through the serialIx of refTree
        save_nodes(Bra, refTree, nodesBra);
        save_nodes(Ket, refTree, nodesKet);
        mrchem::mpi::barrier(mrchem::mpi::comm_orb); // wait until everything is stored before fetching!
    }

    // 3) make dot product for all the nodes and accumulate into S

    int ibank = 0;
#pragma omp parallel for schedule(dynamic) if (serial)
    for (int n = 0; n < max_n; n++) {
        if (n % mpi::orb_size != mpi::orb_rank) continue;
        int csize;
        std::vector<int> orbVecBra; // identifies which Bra orbitals use this node
        std::vector<int> orbVecKet; // identifies which Ket orbitals use this node
        if (parindexVec_ref[n] < 0)
            csize = sizecoeff;
        else
            csize = sizecoeffW;
        if (serial) {
            int node_ix = indexVec_ref[n];      // SerialIx for this node in the reference tree
            int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
            DoubleMatrix coeffBlockBra(csize, node2orbVecBra[node_ix].size());
            DoubleMatrix coeffBlockKet(csize, node2orbVecKet[node_ix].size());
            if (parindexVec_ref[n] < 0) shift = 0;

            for (int j : node2orbVecBra[node_ix]) { // loop over indices of the orbitals using this node
                int orb_node_ix = orb2nodeBra[j][node_ix];
                for (int k = 0; k < csize; k++) coeffBlockBra(k, orbVecBra.size()) = coeffVecBra[j][orb_node_ix][k + shift];
                orbVecBra.push_back(j);
            }
            for (int j : node2orbVecKet[node_ix]) { // loop over indices of the orbitals using this node
                int orb_node_ix = orb2nodeKet[j][node_ix];
                for (int k = 0; k < csize; k++) coeffBlockKet(k, orbVecKet.size()) = coeffVecKet[j][orb_node_ix][k + shift];
                orbVecKet.push_back(j);
            }

            if (orbVecBra.size() > 0 and orbVecKet.size() > 0) {
                DoubleMatrix S_temp(orbVecBra.size(), orbVecKet.size());
                S_temp.noalias() = coeffBlockBra.transpose() * coeffBlockKet;
                for (int i = 0; i < orbVecBra.size(); i++) {
                    for (int j = 0; j < orbVecKet.size(); j++) {
                        if (Bra[orbVecBra[i] % N].spin() == SPIN::Alpha and Ket[orbVecKet[j] % M].spin() == SPIN::Beta) continue;
                        if (Bra[orbVecBra[i] % N].spin() == SPIN::Beta and Ket[orbVecKet[j] % M].spin() == SPIN::Alpha) continue;
                        // must ensure that threads are not competing
                        double &Srealij = Sreal(orbVecBra[i], orbVecKet[j]);
                        double &Stempij = S_temp(i, j);
#pragma omp atomic
                        Srealij += Stempij;
                    }
                }
            }

        } else {
            DoubleMatrix coeffBlockBra(csize, 2 * N);
            DoubleMatrix coeffBlockKet(csize, 2 * M);
            nodesBra.get_nodeblock(indexVec_ref[n], coeffBlockBra.data(), orbVecBra); // get Bra parts
            nodesKet.get_nodeblock(indexVec_ref[n], coeffBlockKet.data(), orbVecKet); // get Ket parts

            if (orbVecBra.size() > 0 and orbVecKet.size() > 0) {
                DoubleMatrix S_temp(orbVecBra.size(), orbVecKet.size());
                coeffBlockBra.conservativeResize(Eigen::NoChange, orbVecBra.size());
                coeffBlockKet.conservativeResize(Eigen::NoChange, orbVecKet.size());
                S_temp.noalias() = coeffBlockBra.transpose() * coeffBlockKet;
                for (int i = 0; i < orbVecBra.size(); i++) {
                    for (int j = 0; j < orbVecKet.size(); j++) {
                        if (Bra[orbVecBra[i] % N].spin() == SPIN::Alpha and Ket[orbVecKet[j] % M].spin() == SPIN::Beta) continue;
                        if (Bra[orbVecBra[i] % N].spin() == SPIN::Beta and Ket[orbVecKet[j] % M].spin() == SPIN::Alpha) continue;
                        Sreal(orbVecBra[i], orbVecKet[j]) += S_temp(i, j);
                    }
                }
            }
        }
    }

    IntVector conjMatBra = IntVector::Zero(N);
    for (int i = 0; i < N; i++) {
        if (!mpi::my_orb(Bra[i])) continue;
        conjMatBra[i] = (Bra[i].conjugate()) ? -1 : 1;
    }
    mpi::allreduce_vector(conjMatBra, mpi::comm_orb);
    IntVector conjMatKet = IntVector::Zero(M);
    for (int i = 0; i < M; i++) {
        if (!mpi::my_orb(Ket[i])) continue;
        conjMatKet[i] = (Ket[i].conjugate()) ? -1 : 1;
    }
    mpi::allreduce_vector(conjMatKet, mpi::comm_orb);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            S.real()(i, j) = Sreal(i, j) + conjMatBra[i] * conjMatKet[j] * Sreal(i + N, j + M);
            S.imag()(i, j) = conjMatKet[j] * Sreal(i, j + M) - conjMatBra[i] * Sreal(i + N, j);
        }
    }

    // 4) collect results from all MPI. Linearity: result is sum of all node contributions

    mpi::allreduce_matrix(S, mpi::comm_orb);

    return S;
}

/** @brief Compute the overlap matrix of the absolute value of the functions S_ij = <|bra_i|||ket_j|>
 *
 */
DoubleMatrix orbital::calc_norm_overlap_matrix(OrbitalVector &BraKet) {
    int N = BraKet.size();
    DoubleMatrix S = DoubleMatrix::Zero(N, N);
    DoubleMatrix Sreal = DoubleMatrix::Zero(2 * N, 2 * N); // same as S, but stored as 4 blocks, rr,ri,ir,ii

    // 1) make union tree without coefficients
    mrcpp::FunctionTree<3> refTree(*MRA);
    mpi::allreduce_Tree_noCoeff(refTree, BraKet, mpi::comm_orb);

    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();

    // get a list of all nodes in union grid, as defined by their indices
    std::vector<double> scalefac;
    std::vector<double *> coeffVec_ref;
    std::vector<int> indexVec_ref;    // serialIx of the nodes
    std::vector<int> parindexVec_ref; // serialIx of the parent nodes
    int max_ix;                       // largest index value (not used here)

    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac, max_ix, refTree);
    int max_n = indexVec_ref.size();

    // only used for serial case:
    std::vector<std::vector<double *>> coeffVec(2 * N);
    std::map<int, std::vector<int>> node2orbVec;     // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2node(2 * N); // for a given orbital and a given node, gives the node index in
                                                     // the orbital given the node index in the reference tree

    bool serial = mpi::orb_size == 1; // flag for serial/MPI switch
    BankAccount nodesBraKet;

    // In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank
    if (serial) {
        // 2) make list of all coefficients, and their reference indices
        // for different orbitals, indexVec will give the same index for the same node in space
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<int> indexVec;    // serialIx of the nodes
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            if (BraKet[j].hasReal()) {
                BraKet[j].real().makeCoeffVector(coeffVec[j], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2node[j][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j);
                }
            }
            if (BraKet[j].hasImag()) {
                BraKet[j].imag().makeCoeffVector(coeffVec[j + N], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2node[j + N][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j + N);
                }
            }
        }
    } else { // MPI case
        // 2) send own nodes to bank, identifying them through the serialIx of refTree
        save_nodes(BraKet, refTree, nodesBraKet);
        mrchem::mpi::barrier(mrchem::mpi::comm_orb); // wait until everything is stored before fetching!
    }

    // 3) make dot product for all the nodes and accumulate into S

    int ibank = 0;
#pragma omp parallel for schedule(dynamic) if (serial)
    for (int n = 0; n < max_n; n++) {
        if (n % mpi::orb_size != mpi::orb_rank) continue;
        int csize;
        int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree
        std::vector<int> orbVec;       // identifies which orbitals use this node
        if (serial and node2orbVec[node_ix].size() <= 0) continue;
        if (parindexVec_ref[n] < 0)
            csize = sizecoeff;
        else
            csize = sizecoeffW;
        // In the serial case we copy the coeff coeffBlock. In the mpi case coeffBlock is provided by the bank
        if (serial) {
            int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
            if (parindexVec_ref[n] < 0) shift = 0;
            DoubleMatrix coeffBlock(csize, node2orbVec[node_ix].size());
            for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                int orb_node_ix = orb2node[j][node_ix];
                for (int k = 0; k < csize; k++) coeffBlock(k, orbVec.size()) = coeffVec[j][orb_node_ix][k + shift];
                orbVec.push_back(j);
            }
            if (orbVec.size() > 0) {
                DoubleMatrix S_temp(orbVec.size(), orbVec.size());
                coeffBlock = coeffBlock.cwiseAbs();
                S_temp.noalias() = coeffBlock.transpose() * coeffBlock;
                for (int i = 0; i < orbVec.size(); i++) {
                    for (int j = 0; j < orbVec.size(); j++) {
                        if (BraKet[orbVec[i] % N].spin() == SPIN::Alpha and BraKet[orbVec[j] % N].spin() == SPIN::Beta) continue;
                        if (BraKet[orbVec[i] % N].spin() == SPIN::Beta and BraKet[orbVec[j] % N].spin() == SPIN::Alpha) continue;
                        double &Srealij = Sreal(orbVec[i], orbVec[j]);
                        double &Stempij = S_temp(i, j);
#pragma omp atomic
                        Srealij += Stempij;
                    }
                }
            }
        } else { // MPI case
            DoubleMatrix coeffBlock(csize, 2 * N);
            nodesBraKet.get_nodeblock(indexVec_ref[n], coeffBlock.data(), orbVec);

            if (orbVec.size() > 0) {
                DoubleMatrix S_temp(orbVec.size(), orbVec.size());
                coeffBlock.conservativeResize(Eigen::NoChange, orbVec.size());
                coeffBlock = coeffBlock.cwiseAbs();
                S_temp.noalias() = coeffBlock.transpose() * coeffBlock;
                for (int i = 0; i < orbVec.size(); i++) {
                    for (int j = 0; j < orbVec.size(); j++) {
                        if (BraKet[orbVec[i] % N].spin() == SPIN::Alpha and BraKet[orbVec[j] % N].spin() == SPIN::Beta) continue;
                        if (BraKet[orbVec[i] % N].spin() == SPIN::Beta and BraKet[orbVec[j] % N].spin() == SPIN::Alpha) continue;
                        Sreal(orbVec[i], orbVec[j]) += S_temp(i, j);
                    }
                }
            }
        }
    }

    IntVector conjMat = IntVector::Zero(N);
    for (int i = 0; i < N; i++) {
        if (!mpi::my_orb(BraKet[i])) continue;
        conjMat[i] = (BraKet[i].conjugate()) ? -1 : 1;
    }
    mpi::allreduce_vector(conjMat, mpi::comm_orb);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            S(i, j) = Sreal(i, j) + conjMat[i] * conjMat[j] * Sreal(i + N, j + N) + conjMat[j] * Sreal(i, j + N) - conjMat[i] * Sreal(i + N, j);
            S(j, i) = S(i, j);
        }
    }

    // Assumes linearity: result is sum of all nodes contributions
    mpi::allreduce_matrix(S, mpi::comm_orb);

    return S;
}

/** @brief Compute Löwdin orthonormalization matrix
 *
 * @param Phi: orbitals to orthonomalize
 *
 * Computes the inverse square root of the orbital overlap matrix S^(-1/2)
 */
ComplexMatrix orbital::calc_lowdin_matrix(OrbitalVector &Phi) {
    Timer overlap_t;
    ComplexMatrix S_tilde = orbital::calc_overlap_matrix(Phi);
    mrcpp::print::time(2, "Computing overlap matrix", overlap_t);
    ComplexMatrix S_m12 = math_utils::hermitian_matrix_pow(S_tilde, -1.0 / 2.0);
    Timer lowdin_t;
    mrcpp::print::time(2, "Computing Lowdin matrix", lowdin_t);
    return S_m12;
}

ComplexMatrix orbital::localize(double prec, OrbitalVector &Phi, ComplexMatrix &F) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Localizing orbitals");
    if (not orbital_vector_is_sane(Phi)) {
        orbital::print(Phi);
        MSG_ABORT("Orbital vector is not sane");
    }
    int nO = Phi.size();
    int nP = size_paired(Phi);
    int nA = size_alpha(Phi);
    int nB = size_beta(Phi);
    ComplexMatrix U = ComplexMatrix::Identity(nO, nO);
    if (nP > 0) U.block(0, 0, nP, nP) = localize(prec, Phi, SPIN::Paired);
    if (nA > 0) U.block(nP, nP, nA, nA) = localize(prec, Phi, SPIN::Alpha);
    if (nB > 0) U.block(nP + nA, nP + nA, nB, nB) = localize(prec, Phi, SPIN::Beta);

    // Transform Fock matrix
    F = U.adjoint() * F * U;
    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Localizing orbitals", t_tot);

    return U;
}

/** @brief Localize a set of orbitals with the same spin

@param Phi_s: Orbital vector containig orbitals with given spin (p/a/b)

Localization is done for each set of spins separately (we don't want to mix spins when localizing).
The localization matrix is returned for further processing.

*/
ComplexMatrix orbital::localize(double prec, OrbitalVector &Phi, int spin) {
    OrbitalVector Phi_s = orbital::disjoin(Phi, spin);
    ComplexMatrix U = calc_localization_matrix(prec, Phi_s);
    Timer rot_t;
    Phi_s = orbital::rotate(Phi_s, U, prec);
    Phi = orbital::adjoin(Phi, Phi_s);
    mrcpp::print::time(2, "Rotating orbitals", rot_t);
    return U;
}

/** @brief Minimize the spatial extension of orbitals, by orbital rotation
 *
 * @param Phi: orbitals to localize (they should all be of the same spin)
 *
 * Minimizes \f$  \sum_{i=1,N}\langle i| {\bf R^2}  | i \rangle - \langle i| {\bf R}| i \rangle^2 \f$
 * which is equivalent to maximizing \f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 *
 * The resulting transformation includes the orthonormalization of the orbitals.
 * Orbitals are rotated in place, and the transformation matrix is returned.
 */
ComplexMatrix orbital::calc_localization_matrix(double prec, OrbitalVector &Phi) {
    ComplexMatrix U;
    int n_it = 0;
    if (Phi.size() > 1) {
        Timer rmat_t;
        RRMaximizer rr(prec, Phi);
        mrcpp::print::time(2, "Computing position matrices", rmat_t);

        Timer rr_t;
        n_it = rr.maximize();
        mrcpp::print::time(2, "Computing Foster-Boys matrix", rr_t);

        if (n_it > 0) {
            println(2, " Foster-Boys localization converged in " << n_it << " iterations!");
            U = rr.getTotalU().cast<ComplexDouble>();
        } else {
            println(2, " Foster-Boys localization did not converge!");
            U = rr.getTotalU().cast<ComplexDouble>();
        }
    } else {
        println(2, " Cannot localize less than two orbitals");
    }
    if (n_it == 0) U = orbital::calc_lowdin_matrix(Phi);
    return U;
}

/** @brief Perform the orbital rotation that diagonalizes the Fock matrix
 *
 * @param Phi: orbitals to rotate
 * @param F: Fock matrix to diagonalize
 *
 * The resulting transformation includes the orthonormalization of the orbitals.
 * Orbitals are rotated in place and Fock matrix is diagonalized in place.
 * The transformation matrix is returned.
 */
ComplexMatrix orbital::diagonalize(double prec, OrbitalVector &Phi, ComplexMatrix &F) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Digonalizing Fock matrix");

    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);
    F = S_m12.adjoint() * F * S_m12;

    Timer diag_t;
    ComplexMatrix U = ComplexMatrix::Zero(F.rows(), F.cols());
    int np = orbital::size_paired(Phi);
    int na = orbital::size_alpha(Phi);
    int nb = orbital::size_beta(Phi);
    if (np > 0) math_utils::diagonalize_block(F, U, 0, np);
    if (na > 0) math_utils::diagonalize_block(F, U, np, na);
    if (nb > 0) math_utils::diagonalize_block(F, U, np + na, nb);
    U = S_m12 * U;
    mrcpp::print::time(2, "Diagonalizing matrix", diag_t);

    Timer rot_t;
    Phi = orbital::rotate(Phi, U, prec);
    mrcpp::print::time(2, "Rotating orbitals", rot_t);

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Diagonalizing Fock matrix", t_tot);
    return U;
}

/** @brief Perform the Löwdin orthonormalization
 *
 * @param Phi: orbitals to orthonormalize
 *
 * Orthonormalizes the orbitals by multiplication of the Löwdin matrix S^(-1/2).
 * Orbitals are rotated in place, and the transformation matrix is returned.
 */
ComplexMatrix orbital::orthonormalize(double prec, OrbitalVector &Phi, ComplexMatrix &F) {
    Timer t_tot, t_lap;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Lowdin orthonormalization");

    ComplexMatrix U = orbital::calc_lowdin_matrix(Phi);

    t_lap.start();
    Phi = orbital::rotate(Phi, U, prec);
    mrcpp::print::time(2, "Rotating orbitals", t_lap);

    // Transform Fock matrix
    F = U.adjoint() * F * U;
    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Lowdin orthonormalization", t_tot);

    return U;
}

/** @brief Returns the number of occupied orbitals */
int orbital::size_occupied(const OrbitalVector &Phi) {
    int nOcc = 0;
    for (auto &phi_i : Phi)
        if (phi_i.occ() > 0) nOcc++;
    return nOcc;
}

/** @brief Returns the number of empty orbitals */
int orbital::size_empty(const OrbitalVector &Phi) {
    int nEmpty = 0;
    for (auto &phi_i : Phi)
        if (phi_i.occ() == 0) nEmpty++;
    return nEmpty;
}

/** @brief Returns the number of singly occupied orbitals */
int orbital::size_singly(const OrbitalVector &Phi) {
    int nSingly = 0;
    for (auto &phi_i : Phi)
        if (phi_i.occ() == 1) nSingly++;
    return nSingly;
}

/** @brief Returns the number of doubly occupied orbitals */
int orbital::size_doubly(const OrbitalVector &Phi) {
    int nDoubly = 0;
    for (auto &phi_i : Phi)
        if (phi_i.occ() == 2) nDoubly++;
    return nDoubly;
}

/** @brief Returns the number of paired orbitals */
int orbital::size_paired(const OrbitalVector &Phi) {
    int nPaired = 0;
    for (auto &phi_i : Phi)
        if (phi_i.spin() == SPIN::Paired) nPaired++;
    return nPaired;
}

/** @brief Returns the number of alpha orbitals */
int orbital::size_alpha(const OrbitalVector &Phi) {
    int nAlpha = 0;
    for (auto &phi_i : Phi)
        if (phi_i.spin() == SPIN::Alpha) nAlpha++;
    return nAlpha;
}

/** @brief Returns the number of beta orbitals */
int orbital::size_beta(const OrbitalVector &Phi) {
    int nBeta = 0;
    for (auto &phi_i : Phi)
        if (phi_i.spin() == SPIN::Beta) nBeta++;
    return nBeta;
}

/** @brief Returns the spin multiplicity of the vector */
int orbital::get_multiplicity(const OrbitalVector &Phi) {
    int nAlpha = get_electron_number(Phi, SPIN::Alpha);
    int nBeta = get_electron_number(Phi, SPIN::Beta);
    int S = std::abs(nAlpha - nBeta);
    return S + 1;
}

/** @brief Returns the number of electrons with the given spin
 *
 * Paired spin (default input) returns the total number of electrons.
 *
 */
int orbital::get_electron_number(const OrbitalVector &Phi, int spin) {
    int nElectrons = 0;
    for (auto &phi_i : Phi) {
        if (spin == SPIN::Paired) {
            nElectrons += phi_i.occ();
        } else if (spin == SPIN::Alpha) {
            if (phi_i.spin() == SPIN::Paired or phi_i.spin() == SPIN::Alpha) nElectrons += 1;
        } else if (spin == SPIN::Beta) {
            if (phi_i.spin() == SPIN::Paired or phi_i.spin() == SPIN::Beta) nElectrons += 1;
        } else {
            MSG_ERROR("Invalid spin argument");
        }
    }
    return nElectrons;
}

/** @brief Returns the total number of nodes in the vector, toggle to get average. */
int orbital::get_n_nodes(const OrbitalVector &Phi, bool avg) {
    long long totNodes = 0;
    for (const auto &phi_i : Phi) totNodes += phi_i.getNNodes(NUMBER::Total);
    if (avg) totNodes /= Phi.size();
    if (totNodes > INT_MAX) MSG_WARN("Integer overflow: " << totNodes);
    return static_cast<int>(totNodes);
}

/** @brief Returns the size of the coefficients of all nodes in the vector in kBytes, toggle to get average.*/
int orbital::get_size_nodes(const OrbitalVector &Phi, bool avg) {
    long long totSize = 0;
    for (const auto &phi_i : Phi) totSize += phi_i.getSizeNodes(NUMBER::Total);
    if (avg) totSize /= Phi.size();
    if (totSize > INT_MAX) MSG_WARN("Integer overflow: " << totSize);
    return static_cast<int>(totSize);
}

/** @brief Returns a vector containing the orbital spins */
IntVector orbital::get_spins(const OrbitalVector &Phi) {
    int nOrbs = Phi.size();
    IntVector spins = IntVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) spins(i) = Phi[i].spin();
    return spins;
}

/** @brief Assigns spin to each orbital
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void orbital::set_spins(OrbitalVector &Phi, const IntVector &spins) {
    if (Phi.size() != spins.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < Phi.size(); i++) Phi[i].setSpin(spins(i));
}

/** @brief Returns a vector containing the orbital occupations */
IntVector orbital::get_occupations(const OrbitalVector &Phi) {
    int nOrbs = Phi.size();
    IntVector occ = IntVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) occ(i) = Phi[i].occ();
    return occ;
}

/** @brief Assigns occupation to each orbital
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void orbital::set_occupations(OrbitalVector &Phi, const IntVector &occ) {
    if (Phi.size() != occ.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < Phi.size(); i++) Phi[i].setOcc(occ(i));
}

/** @brief Returns a vector containing the orbital square norms */
DoubleVector orbital::get_squared_norms(const OrbitalVector &Phi) {
    int nOrbs = Phi.size();
    DoubleVector norms = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        if (mpi::my_orb(Phi[i])) norms(i) = Phi[i].squaredNorm();
    }
    mpi::allreduce_vector(norms, mpi::comm_orb);
    return norms;
}

/** @brief Returns a vector containing the orbital norms */
DoubleVector orbital::get_norms(const OrbitalVector &Phi) {
    int nOrbs = Phi.size();
    DoubleVector norms = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        if (mpi::my_orb(Phi[i])) norms(i) = Phi[i].norm();
    }
    mpi::allreduce_vector(norms, mpi::comm_orb);
    return norms;
}

/** @brief Returns a vector containing the orbital integrals */
ComplexVector orbital::get_integrals(const OrbitalVector &Phi) {
    int nOrbs = Phi.size();
    ComplexVector ints = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        if (mpi::my_orb(Phi[i])) ints(i) = Phi[i].integrate();
    }
    mpi::allreduce_vector(ints, mpi::comm_orb);
    return ints;
}

/** @brief Checks if a vector of orbitals is correctly ordered (paired/alpha/beta) */
bool orbital::orbital_vector_is_sane(const OrbitalVector &Phi) {
    int nO = Phi.size();
    int nP = size_paired(Phi);
    int nA = size_alpha(Phi);
    int nB = size_beta(Phi);
    int previous_spin = 0;

    if (nO != nP + nA + nB) return false; // not all orbitals are accounted for

    for (int i = 0; i < nO; i++) {
        if (Phi[i].spin() < previous_spin) return false; // wrong orbital order
        previous_spin = Phi[i].spin();
    }
    return true; // sane orbital set
}
/** @brief Returns the start index of a given orbital type (p/a/b)
 *
 *  Returns a negative number if the type of orbitals is not present.
 *  The ordering of orbitals in a given OrbitalVector is fixed and
 *  this can be used to determine the end index as well.
 */
int orbital::start_index(const OrbitalVector &Phi, int spin) {
    int nOrbs = Phi.size();
    for (int i = 0; i < nOrbs; i++) {
        if (Phi[i].spin() == spin) return i;
    }
    return -1;
}

void orbital::print(const OrbitalVector &Phi) {
    auto pprec = Printer::getPrecision();
    auto w0 = Printer::getWidth() - 1;
    auto w1 = 5;
    auto w2 = 2 * w0 / 9;
    auto w3 = w0 - 3 * w1 - 3 * w2;

    auto N_e = orbital::get_electron_number(Phi);
    auto N_a = orbital::size_alpha(Phi) + orbital::size_paired(Phi);
    auto N_b = orbital::size_beta(Phi) + orbital::size_paired(Phi);

    std::stringstream o_head;
    o_head << std::setw(w1) << "n";
    o_head << std::setw(w1) << "Occ";
    o_head << std::setw(w1) << "Spin";
    o_head << std::string(w3 - 1, ' ') << ':';
    o_head << std::setw(3 * w2) << "Norm";

    mrcpp::print::header(0, "Molecular Orbitals");
    print_utils::scalar(0, "Alpha electrons ", N_a, "", 0, false);
    print_utils::scalar(0, "Beta electrons  ", N_b, "", 0, false);
    print_utils::scalar(0, "Total electrons ", N_e, "", 0, false);
    mrcpp::print::separator(0, '-');
    println(0, o_head.str());
    mrcpp::print::separator(0, '-');

    auto norms = orbital::get_norms(Phi); // includes allreduce

    auto nodes = 0;
    auto memory = 0.0;
    for (int i = 0; i < Phi.size(); i++) {
        nodes += Phi[i].getNNodes(NUMBER::Total);
        memory += Phi[i].getSizeNodes(NUMBER::Total) / 1024.0;
        std::stringstream o_txt;
        o_txt << std::setw(w1 - 1) << i;
        o_txt << std::setw(w1) << Phi[i].occ();
        o_txt << std::setw(w1) << Phi[i].printSpin();
        print_utils::scalar(0, o_txt.str(), norms[i], "", 2 * pprec, true);
    }

    mrcpp::print::separator(2, '-');
    print_utils::scalar(2, "Total MO nodes ", nodes, "", 0, false);
    print_utils::scalar(2, "Total MO memory ", memory, "(MB)", 2, false);
    mrcpp::print::separator(0, '=', 2);
}

DoubleVector orbital::calc_eigenvalues(const OrbitalVector &Phi, const ComplexMatrix &F_mat) {
    if (F_mat.cols() != Phi.size()) MSG_ABORT("Invalid Fock matrix");
    if (not orbital::orbital_vector_is_sane(Phi)) MSG_ABORT("Insane orbital vector");

    DoubleVector epsilon = DoubleVector::Zero(Phi.size());
    int np = orbital::size_paired(Phi);
    int na = orbital::size_alpha(Phi);
    int nb = orbital::size_beta(Phi);
    if (np > 0) {
        Timer timer;
        Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(np);
        es.compute(F_mat.block(0, 0, np, np));
        epsilon.segment(0, np) = es.eigenvalues();
        mrcpp::print::time(1, "Diagonalize Fock matrix", timer);
    }
    if (na > 0) {
        Timer timer;
        Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(na);
        es.compute(F_mat.block(np, np, na, na));
        epsilon.segment(np, na) = es.eigenvalues();
        mrcpp::print::time(1, "Diagonalize Fock matrix (alpha)", timer);
    }
    if (nb > 0) {
        Timer timer;
        Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(nb);
        es.compute(F_mat.block(np + na, np + na, nb, nb));
        epsilon.segment(np + na, nb) = es.eigenvalues();
        mrcpp::print::time(1, "Diagonalize Fock matrix (beta)", timer);
    }
    return epsilon;
}

/** @brief Prints statistics about the size of orbitals in an OrbitalVector
 *
 * This is a collective function. Can be made non-collective by setting all = false.
 * outputs respectively:
 * Total size of orbital vector, average per MPI, Max per MPI, Max (largest)
 * orbital, smallest orbital, max total (not only the orbitalvector) memory
 * usage among all MP, minimum total (not only the orbitalvector) memory
 * usage among all MPI
 *
 */
int orbital::print_size_nodes(const OrbitalVector &Phi, const std::string &txt, bool all, int plevel) {
    double nMax = 0.0, vMax = 0.0; // node max, vector max
    double nMin = 9.9e9, vMin = 9.9e9;
    double nSum = 0.0, vSum = 0.0;
    double nOwnOrbs = 0.0, ownSumMax = 0.0, ownSumMin = 9.9e9;
    double totMax = 0.0, totMin = 9.9e9;
    println(0, "OrbitalVector sizes statistics " << txt << " (MB)");

    IntVector sNodes = IntVector::Zero(Phi.size());
    for (int i = 0; i < Phi.size(); i++) sNodes[i] = Phi[i].getSizeNodes(NUMBER::Total);

    // stats for own orbitals
    for (int i = 0; i < Phi.size(); i++) {
        if (sNodes[i] > 0) {
            nOwnOrbs++;
            if (sNodes[i] > nMax) nMax = sNodes[i];
            if (sNodes[i] < nMin) nMin = sNodes[i];
            nSum += sNodes[i];
        }
    }
    if (nSum == 0.0) nMin = 0.0;

    DoubleMatrix vecStats = DoubleMatrix::Zero(5, mpi::orb_size);
    vecStats(0, mpi::orb_rank) = nMax;
    vecStats(1, mpi::orb_rank) = nMin;
    vecStats(2, mpi::orb_rank) = nSum;
    vecStats(3, mpi::orb_rank) = nOwnOrbs;
    vecStats(4, mpi::orb_rank) = mrcpp::details::get_memory_usage();

    if (all) {
        mpi::allreduce_matrix(vecStats, mpi::comm_orb);
        // overall stats
        for (int i = 0; i < mpi::orb_size; i++) {
            if (vecStats(0, i) > vMax) vMax = vecStats(0, i);
            if (vecStats(1, i) < vMin) vMin = vecStats(1, i);
            if (vecStats(2, i) > ownSumMax) ownSumMax = vecStats(2, i);
            if (vecStats(2, i) < ownSumMin) ownSumMin = vecStats(2, i);
            if (vecStats(4, i) > totMax) totMax = vecStats(4, i);
            if (vecStats(4, i) < totMin) totMin = vecStats(4, i);
            vSum += vecStats(2, i);
        }
    } else {
        int i = mpi::orb_rank;
        if (vecStats(0, i) > vMax) vMax = vecStats(0, i);
        if (vecStats(1, i) < vMin) vMin = vecStats(1, i);
        if (vecStats(2, i) > ownSumMax) ownSumMax = vecStats(2, i);
        if (vecStats(2, i) < ownSumMin) ownSumMin = vecStats(2, i);
        if (vecStats(4, i) > totMax) totMax = vecStats(4, i);
        if (vecStats(4, i) < totMin) totMin = vecStats(4, i);
        vSum += vecStats(2, i);
    }
    totMax /= 1024.0;
    totMin /= 1024.0;
    printout(plevel, "Total orbvec " << static_cast<int>(vSum / 1024));
    printout(plevel, ", Av/MPI " << static_cast<int>(vSum / 1024 / mpi::orb_size));
    printout(plevel, ", Max/MPI " << static_cast<int>(ownSumMax / 1024));
    printout(plevel, ", Max/orb " << static_cast<int>(vMax / 1024));
    printout(plevel, ", Min/orb " << static_cast<int>(vMin / 1024));

    auto totMinInt = static_cast<int>(totMin);
    auto totMaxInt = static_cast<int>(totMax);
    if (all) {
        println(plevel, ", Total max " << totMaxInt << ", Total min " << totMinInt << " MB");
    } else {
        println(plevel, ", Total master " << totMaxInt << " MB");
    }
    return vSum;
}

} // namespace mrchem
