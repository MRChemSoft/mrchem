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

#include <fstream>

#include <MRCPP/Printer>
#include <MRCPP/Timer>
#include <MRCPP/trees/FunctionNode.h>
#include <MRCPP/utils/details.h>

#include "utils/RRMaximizer.h"
#include "utils/math_utils.h"
#include "utils/print_utils.h"

#include "Orbital.h"
#include "orbital_utils.h"

using mrcpp::FunctionNode;
using mrcpp::FunctionTreeVector;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

namespace orbital {
ComplexMatrix localize(double prec, OrbitalVector &Phi, int spin);
ComplexMatrix calc_localization_matrix(double prec, OrbitalVector &Phi);

/* POD struct for orbital meta data. Used for simple MPI communication. */
struct OrbitalData {
    int rank_id;
    int spin;
    double occ;
};
OrbitalData getOrbitalData(const Orbital &orb) {
    OrbitalData orb_data;
    orb_data.rank_id = orb.getRank();
    orb_data.spin = orb.spin();
    orb_data.occ = orb.occ();
    return orb_data;
}
} // namespace orbital

/****************************************
 * Orbital related standalone functions *
 ****************************************/

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
        if (mrcpp::mpi::my_func(Phi_a[i]) != mrcpp::mpi::my_func(Phi_b[i])) MSG_ABORT("MPI rank mismatch");
        mrcpp::add(out[i], a, Phi_a[i], b, Phi_b[i], prec);
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
    OrbitalVector Psi = orbital::deep_copy(Phi);
    mrcpp::rotate(Psi, U, prec);
    return Psi;
}

/** @brief Save all nodes in bank; identify them using serialIx from refTree
 * shift is a shift applied in the id
 */
void orbital::save_nodes(OrbitalVector Phi, mrcpp::FunctionTree<3> &refTree, mrcpp::BankAccount &account) {
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
        if (not mrcpp::mpi::my_func(Phi[j])) continue;
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
        Orbital out_i;
        if (mrcpp::mpi::my_func(i)) mrcpp::deep_copy(out_i, i);
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
        Orbital out_i;
        out_i.func_ptr->data = i.func_ptr->data;
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
    for (auto &phi : Phi_a) {
        if (phi.getRank() % mrcpp::mpi::wrk_size != out.size() % mrcpp::mpi::wrk_size) {
            // need to send orbital from owner to new owner
            if (mrcpp::mpi::my_func(phi)) { mrcpp::mpi::send_function(phi, out.size() % mrcpp::mpi::wrk_size, phi.getRank(), mrcpp::mpi::comm_wrk); }
            if (mrcpp::mpi::my_func(out.size())) { mrcpp::mpi::recv_function(phi, phi.getRank() % mrcpp::mpi::wrk_size, phi.getRank(), mrcpp::mpi::comm_wrk); }
        }
        phi.setRank(out.size());
        out.push_back(phi);
    }
    for (auto &phi : Phi_b) {
        if (phi.getRank() % mrcpp::mpi::wrk_size != out.size() % mrcpp::mpi::wrk_size) {
            // need to send orbital from owner to new owner
            if (mrcpp::mpi::my_func(phi)) { mrcpp::mpi::send_function(phi, out.size() % mrcpp::mpi::wrk_size, phi.getRank(), mrcpp::mpi::comm_wrk); }
            if (mrcpp::mpi::my_func(out.size())) { mrcpp::mpi::recv_function(phi, phi.getRank() % mrcpp::mpi::wrk_size, phi.getRank(), mrcpp::mpi::comm_wrk); }
        }
        phi.setRank(out.size());
        out.push_back(phi);
    }
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
    for (auto &Phi_i : Phi) {
        Orbital i(Phi_i);
        if (i.spin() == spin) {
            if (i.getRank() % mrcpp::mpi::wrk_size != out.size() % mrcpp::mpi::wrk_size) {
                // need to send orbital from owner to new owner
                if (mrcpp::mpi::my_func(i)) { mrcpp::mpi::send_function(i, out.size() % mrcpp::mpi::wrk_size, i.getRank(), mrcpp::mpi::comm_wrk); }
                if (mrcpp::mpi::my_func(out.size())) { mrcpp::mpi::recv_function(i, i.getRank() % mrcpp::mpi::wrk_size, i.getRank(), mrcpp::mpi::comm_wrk); }
            }
            i.setRank(out.size());
            out.push_back(i);
        } else {
            if (i.getRank() % mrcpp::mpi::wrk_size != tmp.size() % mrcpp::mpi::wrk_size) {
                // need to send orbital from owner to new owner
                if (mrcpp::mpi::my_func(i)) { mrcpp::mpi::send_function(i, tmp.size() % mrcpp::mpi::wrk_size, i.getRank(), mrcpp::mpi::comm_wrk); }
                if (mrcpp::mpi::my_func(tmp.size())) { mrcpp::mpi::recv_function(i, i.getRank() % mrcpp::mpi::wrk_size, i.getRank(), mrcpp::mpi::comm_wrk); }
            }
            i.setRank(tmp.size());
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
void orbital::save_orbitals(OrbitalVector &Phi, const std::string &file, int spin, int text_format) {
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
            if (mrcpp::mpi::my_func(Phi[i])) saveOrbital(orbname.str(), Phi[i], text_format);
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
        loadOrbital(orbname.str(), phi_i);
        phi_i.setRank(i);
        if (phi_i.hasReal() or phi_i.hasImag()) {
            Phi.push_back(phi_i);
            print_utils::qmfunction(2, "'" + orbname.str() + "'", phi_i, t1);
            if (not mrcpp::mpi::my_func(phi_i)) phi_i.free();
        } else {
            break;
        }
    }
    mrcpp::print::footer(2, t_tot, 2);
    return Phi;
}

/** @brief Normalize single orbital. Private function. */
void orbital::normalize(Orbital phi) {
    phi.rescale(1.0 / phi.norm());
}

/** @brief Normalize all orbitals in the set */
void orbital::normalize(OrbitalVector &Phi) {
    mrcpp::mpi::free_foreign(Phi);
    for (auto &phi_i : Phi)
        if (mrcpp::mpi::my_func(phi_i)) orbital::normalize(phi_i);
}

/** @brief In place orthogonalize against inp. Private function. */
void orbital::orthogonalize(double prec, Orbital &&phi, Orbital psi) {
    ComplexDouble overlap = mrcpp::dot(psi, phi);
    double sq_norm = psi.getSquareNorm();
    if (std::abs(overlap) > prec) phi.add(-1.0 * overlap / sq_norm, psi);
}

/** @brief Gram-Schmidt orthogonalize orbitals within the set */
void orbital::orthogonalize(double prec, OrbitalVector &Phi) {
    mrcpp::mpi::free_foreign(Phi);
    for (int i = 0; i < Phi.size(); i++) {
        for (int j = 0; j < i; j++) {
            int tag = 7632 * i + j;
            int src = (Phi[j].getRank()) % mrcpp::mpi::wrk_size;
            int dst = (Phi[i].getRank()) % mrcpp::mpi::wrk_size;
            if (mrcpp::mpi::my_func(Phi[i]) and mrcpp::mpi::my_func(Phi[j])) {
                mrcpp::orthogonalize(prec / Phi.size(), Phi[i], Phi[j]);
            } else {
                if (mrcpp::mpi::my_func(Phi[i])) {
                    mrcpp::mpi::recv_function(Phi[j], src, tag, mrcpp::mpi::comm_wrk);
                    mrcpp::orthogonalize(prec / Phi.size(), Phi[i], Phi[j]);
                    Phi[j].free();
                }
                if (mrcpp::mpi::my_func(Phi[j])) mrcpp::mpi::send_function(Phi[j], dst, tag, mrcpp::mpi::comm_wrk);
            }
        }
    }
}

OrbitalChunk orbital::get_my_chunk(OrbitalVector &Phi) {
    OrbitalChunk chunk;
    for (int i = 0; i < Phi.size(); i++) {
        if (mrcpp::mpi::my_func(i)) chunk.push_back(std::make_tuple(i, Orbital(Phi[i])));
    }
    return chunk;
}

/** @brief Orthogonalize the Phi orbitals against all orbitals in Psi.
 *  orthogonal spins means orthogonal orbitals.
 */
void orbital::orthogonalize(double prec, OrbitalVector &Phi, OrbitalVector &Psi) {
    mrcpp::orthogonalize(prec, Phi, Psi);
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
    return mrcpp::calc_overlap_matrix(BraKet);
}

/** @brief Compute the overlap matrix S_ij = <bra_i|ket_j>
 *
 */
ComplexMatrix orbital::calc_overlap_matrix(OrbitalVector &Bra, OrbitalVector &Ket) {
    return mrcpp::calc_overlap_matrix(Bra, Ket);
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
    Timer lowdin_t;
    ComplexMatrix S_m12 = math_utils::hermitian_matrix_pow(S_tilde, -1.0 / 2.0);
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
    mrcpp::rotate(Phi_s, U, prec);
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
    mrcpp::rotate(Phi, U, prec);
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
    mrcpp::rotate(Phi, U, prec);
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
    double nAlpha = get_electron_number(Phi, SPIN::Alpha);
    double nBeta = get_electron_number(Phi, SPIN::Beta);
    // LR: not sure this is currently meaningful for fractional occupancies
    // as a minimum should we add some kind of warning/error message if the result is not an integer?
    int S = std::abs(nAlpha - nBeta);
    return S + 1;
}

/** @brief Returns the number of electrons with the given spin
 *
 * Paired spin (default input) returns the total number of electrons.
 *
 */
double orbital::get_electron_number(const OrbitalVector &Phi, int spin) {
    double nElectrons = 0;
    for (auto &phi_i : Phi) {
        if (spin == SPIN::Paired) {
            nElectrons += (int)phi_i.occ() + 0.5; // nearest integer
        } else if (spin == SPIN::Alpha) {
            if (phi_i.spin() == SPIN::Paired) {nElectrons += 0.5 * phi_i.occ();}
            else if (phi_i.spin() == SPIN::Alpha) {nElectrons += phi_i.occ();}
        } else if (spin == SPIN::Beta) {
            if (phi_i.spin() == SPIN::Paired) {nElectrons += 0.5 * phi_i.occ();}
            else if (phi_i.spin() == SPIN::Beta) {nElectrons += phi_i.occ();}
        } else {
            MSG_ERROR("Invalid spin argument");
        }
    }
    return nElectrons;
}

/** @brief Returns the total number of nodes in the vector, toggle to get average. */
int orbital::get_n_nodes(const OrbitalVector &Phi, bool avg) {
    long long totNodes = 0;
    int mysize = 0;
    for (const auto &phi_i : Phi) totNodes += phi_i.getNNodes();
    for (const auto &phi_i : Phi)
        if (mrcpp::mpi::my_func(phi_i)) mysize++;
    if (avg and mysize > 0) totNodes /= mysize;
    if (totNodes > INT_MAX) MSG_WARN("Integer overflow: " << totNodes);
    return static_cast<int>(totNodes);
}

/** @brief Returns the size of the coefficients of all nodes in the vector in kBytes, toggle to get average.*/
int orbital::get_size_nodes(const OrbitalVector &Phi, bool avg) {
    long long totSize = 0;
    int mysize = 0;
    for (const auto &phi_i : Phi) totSize += phi_i.getSizeNodes();
    for (const auto &phi_i : Phi)
        if (mrcpp::mpi::my_func(phi_i)) mysize++;
    if (avg and mysize > 0) totSize /= mysize;
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
    for (int i = 0; i < Phi.size(); i++) Phi[i].spin() = i;
}

/** @brief Returns a vector containing the orbital occupations */
DoubleVector orbital::get_occupations(const OrbitalVector &Phi) {
    int nOrbs = Phi.size();
    DoubleVector occup = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) occup(i) = Phi[i].occ();
    return occup;
}

/** @brief Assigns occupation to each orbital
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void orbital::set_occupations(OrbitalVector &Phi, const DoubleVector &occup) {
    if (Phi.size() != occup.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < Phi.size(); i++) Phi[i].occ() = occup(i);
}

/** @brief Returns a vector containing the orbital square norms */
DoubleVector orbital::get_squared_norms(const OrbitalVector &Phi) {
    int nOrbs = Phi.size();
    DoubleVector norms = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        if (mrcpp::mpi::my_func(Phi[i])) norms(i) = Phi[i].getSquareNorm();
    }
    mrcpp::mpi::allreduce_vector(norms, mrcpp::mpi::comm_wrk);
    return norms;
}

/** @brief Returns a vector containing the orbital norms */
DoubleVector orbital::get_norms(const OrbitalVector &Phi) {
    int nOrbs = Phi.size();
    DoubleVector norms = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        if (mrcpp::mpi::my_func(Phi[i])) norms(i) = Phi[i].norm();
    }
    mrcpp::mpi::allreduce_vector(norms, mrcpp::mpi::comm_wrk);
    return norms;
}

/** @brief Returns a vector containing the orbital integrals */
ComplexVector orbital::get_integrals(const OrbitalVector &Phi) {
    int nOrbs = Phi.size();
    ComplexVector ints = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        if (mrcpp::mpi::my_func(Phi[i])) ints(i) = Phi[i].integrate();
    }
    mrcpp::mpi::allreduce_vector(ints, mrcpp::mpi::comm_wrk);
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
    auto N_a = orbital::get_electron_number(Phi, SPIN::Alpha);
    auto N_b = orbital::get_electron_number(Phi, SPIN::Beta);

    std::stringstream o_head;
    o_head << std::setw(w1) << "n";
    o_head << std::setw(w1 + 3) << "Occ";
    o_head << std::setw(w1) << "Spin";
    o_head << std::string(w3 - 4, ' ') << ':';
    o_head << std::setw(3 * w2) << "Norm";

    mrcpp::print::header(0, "Molecular Orbitals");
    print_utils::scalar(0, "Alpha electrons ", N_a, "", 3, false);
    print_utils::scalar(0, "Beta electrons  ", N_b, "", 3, false);
    print_utils::scalar(0, "Total electrons ", N_e, "", 3, false);
    mrcpp::print::separator(0, '-');
    println(0, o_head.str());
    mrcpp::print::separator(0, '-');

    auto norms = orbital::get_norms(Phi); // includes allreduce

    auto nodes = 0;
    auto memory = 0.0;
    for (int i = 0; i < Phi.size(); i++) {
        nodes += Phi[i].getNNodes();
        memory += Phi[i].getSizeNodes() / 1024.0;
        std::stringstream o_txt;
        o_txt << std::setw(w1 - 1) << i;
        o_txt << std::setw(w1 + 3) << std::setprecision(2) << std::fixed << Phi[i].occ();
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
    for (int i = 0; i < Phi.size(); i++) sNodes[i] = Phi[i].getSizeNodes();

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

    DoubleMatrix vecStats = DoubleMatrix::Zero(5, mrcpp::mpi::wrk_size);
    vecStats(0, mrcpp::mpi::wrk_rank) = nMax;
    vecStats(1, mrcpp::mpi::wrk_rank) = nMin;
    vecStats(2, mrcpp::mpi::wrk_rank) = nSum;
    vecStats(3, mrcpp::mpi::wrk_rank) = nOwnOrbs;
    vecStats(4, mrcpp::mpi::wrk_rank) = mrcpp::details::get_memory_usage();

    if (all) {
        mrcpp::mpi::allreduce_matrix(vecStats, mrcpp::mpi::comm_wrk);
        // overall stats
        for (int i = 0; i < mrcpp::mpi::wrk_size; i++) {
            if (vecStats(0, i) > vMax) vMax = vecStats(0, i);
            if (vecStats(1, i) < vMin) vMin = vecStats(1, i);
            if (vecStats(2, i) > ownSumMax) ownSumMax = vecStats(2, i);
            if (vecStats(2, i) < ownSumMin) ownSumMin = vecStats(2, i);
            if (vecStats(4, i) > totMax) totMax = vecStats(4, i);
            if (vecStats(4, i) < totMin) totMin = vecStats(4, i);
            vSum += vecStats(2, i);
        }
    } else {
        int i = mrcpp::mpi::wrk_rank;
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
    printout(plevel, ", Av/MPI " << static_cast<int>(vSum / 1024 / mrcpp::mpi::wrk_size));
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

void orbital::saveOrbital(const std::string &file, mrcpp::CompFunction<3> &orb, int text_format) {
    orbital::saveOrbital(file, Orbital(orb), text_format);
}

/** @brief Write orbital to disk
 *
 * @param file: file name prefix
 *
 * Given a file name prefix (e.g. "phi_0"), this will produce separate
 * binary files for meta data ("phi_0.meta"), real ("phi_0_re.tree")
 * and imaginary ("phi_0_im.tree") parts.
 */
void orbital::saveOrbital(const std::string &file, const Orbital &orb, int text_format) {
    // writing meta data
    std::stringstream metafile;
    metafile << file << ".meta";

    if (not text_format) {
        std::fstream f;
        f.open(metafile.str(), std::ios::out | std::ios::binary);
        if (not f.is_open()) MSG_ERROR("Unable to open file");
        mrcpp::CompFunctionData<3> orbdata = orb.getFuncData();
        f.write((char *)&orb.func_ptr->data, sizeof(mrcpp::CompFunctionData<3>));
        f.close();
    }

    // writing real tree
    if (orb.isreal()) {
        std::stringstream fname;
        fname << file << "_real";
        if (text_format)
            orb.CompD[0]->saveTreeTXT(fname.str());
        else
            orb.CompD[0]->saveTree(fname.str());
    }

    // writing complex tree
    if (orb.iscomplex()) {
        std::stringstream fname;
        fname << file << "_complex";
        if (text_format)
            orb.CompC[0]->saveTreeTXT(fname.str());
        else
            orb.CompC[0]->saveTree(fname.str());
    }
}

/** @brief Read orbital from disk
 *
 * @param file: file name prefix
 *
 * Given a file name prefix (e.g. "phi_0"), this will read separate
 * binary files for meta data ("phi_0.meta"), real ("phi_0_re.tree")
 * and imaginary ("phi_0_im.tree") parts.
 */
void orbital::loadOrbital(const std::string &file, Orbital &orb) {
    if (orb.hasReal()) MSG_ERROR("Orbital not empty");
    if (orb.hasImag()) MSG_ERROR("Orbital not empty");

    // first test if the file is in text format or MRChem binary format
    std::ifstream testfile;
    std::stringstream fname_re;
    fname_re << file << "_real";
    testfile.open(fname_re.str());
    if (testfile) {
        // since the MRChem file names end by .tree, we assume that this one is in text format
        orb.defreal();
        orb.alloc(1);
        orb.CompD[0]->loadTreeTXT(fname_re.str());
        return;
    }
    std::stringstream fname_co;
    fname_co << file << "_complex";
    testfile.open(fname_co.str());
    if (testfile) {
        // since the MRChem file names end by .tree, we assume that this one is in text format
        orb.defcomplex();
        orb.alloc(1);
        orb.CompC[0]->loadTreeTXT(fname_co.str());
        return;
    }

    // the file is not in TXT format if this point is reached

    // reading meta data
    std::stringstream fmeta;
    fmeta << file << ".meta";

    std::fstream f;
    f.open(fmeta.str(), std::ios::in | std::ios::binary);
    if (f.is_open()) f.read((char *)&orb.func_ptr->data, sizeof(mrcpp::CompFunctionData<3>));
    f.close();

    std::array<int, 3> corner{orb.data().corner[0], orb.data().corner[1], orb.data().corner[2]};
    std::array<int, 3> boxes{orb.data().boxes[0], orb.data().boxes[1], orb.data().boxes[2]};
    mrcpp::BoundingBox<3> world(orb.data().scale, corner, boxes);

    mrcpp::MultiResolutionAnalysis<3> *mra = nullptr;
    if (orb.data().type == mrcpp::Interpol) {
        mrcpp::InterpolatingBasis basis(orb.data().order);
        mra = new mrcpp::MultiResolutionAnalysis<3>(world, basis, orb.data().depth);
    } else if (orb.data().type == mrcpp::Legendre) {
        mrcpp::LegendreBasis basis(orb.data().order);
        mra = new mrcpp::MultiResolutionAnalysis<3>(world, basis, orb.data().depth);
    } else {
        MSG_ABORT("Invalid basis type!");
    }

    // reading real orbital
    if (orb.isreal()) {
        std::stringstream fname;
        fname << file << "_real";
        orb.alloc(1);
        orb.CompD[0]->loadTree(fname.str());
    }

    // reading complex orbital
    if (orb.iscomplex()) {
        std::stringstream fname;
        fname << file << "_complex";
        orb.alloc(1);
        orb.CompC[0]->loadTree(fname.str());
    }
    delete mra;
}

/** @brief Returns a character representing the spin (a/b/p) */
// char orbital::printSpin(const Orbital& orb) {
//    char sp = 'u';
//    if (orb.spin() == SPIN::Paired) sp = 'p';
//    if (orb.spin() == SPIN::Alpha) sp = 'a';
//    if (orb.spin() == SPIN::Beta) sp = 'b';
//    return sp;
//}

} // namespace mrchem
