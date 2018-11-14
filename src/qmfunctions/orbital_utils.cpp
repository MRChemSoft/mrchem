/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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

#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"
#include "utils/RRMaximizer.h"
#include "utils/math_utils.h"

#include "Orbital.h"
#include "OrbitalIterator.h"
#include "orbital_utils.h"
#include "qmfunction_utils.h"

using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

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
    if (Bra.size() != Ket.size()) MSG_FATAL("Size mismatch");

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
        if (not mpi::my_orb(Bra[i])) Bra[i].free();
    }
    mpi::allreduce_vector(result, mpi::comm_orb);
    return result;
}

/** @brief Compare spin and occupancy of two orbitals
 *
 *  Returns true if orbital parameters are the same.
 *
 */
bool orbital::compare(const Orbital &phi_a, const Orbital &phi_b) {
    bool comp = true;
    if (compare_occ(phi_a, phi_b) < 0) {
        MSG_WARN("Different occupancy");
        comp = false;
    }
    if (compare_spin(phi_a, phi_b) < 0) {
        MSG_WARN("Different spin");
        comp = false;
    }
    return comp;
}

/** @brief Compare occupancy of two orbitals
 *
 *  Returns the common occupancy if they match, -1 if they differ.
 *
 */
int orbital::compare_occ(const Orbital &phi_a, const Orbital &phi_b) {
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

/** @brief out_i = a*(inp_a)_i + b*(inp_b)_i
 *
 *  Component-wise addition of orbitals.
 *
 */
OrbitalVector orbital::add(ComplexDouble a, OrbitalVector &Phi_a, ComplexDouble b, OrbitalVector &Phi_b, double prec) {
    if (Phi_a.size() != Phi_b.size()) MSG_ERROR("Size mismatch");

    OrbitalVector out = orbital::param_copy(Phi_a);
    for (int i = 0; i < Phi_a.size(); i++) {
        if (Phi_a[i].rankID() != Phi_b[i].rankID()) MSG_FATAL("MPI rank mismatch");
        qmfunction::add(out[i], a, Phi_a[i], b, Phi_b[i], prec);
    }
    return out;
}

/** @brief Orbital transformation out_vec = U*inp_vec
 *
 * MPI: Rank distribution of output vector is the same as input vector
 *
 */
OrbitalVector orbital::rotate(const ComplexMatrix &U, OrbitalVector &Phi, double prec) {
    // Get all out orbitals belonging to this MPI
    double inter_prec = (mpi::numerically_exact) ? -1.0 : prec;
    OrbitalVector out = orbital::param_copy(Phi);
    OrbitalIterator iter(Phi);
    while (iter.next()) {
        for (int i = 0; i < out.size(); i++) {
            if (not mpi::my_orb(out[i])) continue;
            ComplexVector coef_vec(iter.get_size());
            QMFunctionVector func_vec;
            for (int j = 0; j < iter.get_size(); j++) {
                int idx_j = iter.idx(j);
                Orbital &recv_j = iter.orbital(j);
                coef_vec[j] = U(i, idx_j);
                func_vec.push_back(recv_j);
            }
            Orbital tmp_i = out[i].paramCopy();
            qmfunction::linear_combination(tmp_i, coef_vec, func_vec, inter_prec);
            out[i].add(1.0, tmp_i); // In place addition
            out[i].crop(inter_prec);
            tmp_i.free();
        }
    }

    if (mpi::numerically_exact) {
        for (auto &out_i : out) {
            if (mpi::my_orb(out_i)) out_i.crop(prec);
        }
    }

    return out;
}

/** @brief Deep copy
 *
 * New orbitals are constructed as deep copies of the input set.
 *
 */
OrbitalVector orbital::deep_copy(OrbitalVector &Phi) {
    OrbitalVector out;
    for (int i = 0; i < Phi.size(); i++) {
        Orbital out_i = Phi[i].deepCopy();
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
    for (int i = 0; i < Phi.size(); i++) {
        Orbital out_i = Phi[i].paramCopy();
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
    for (int i = 0; i < Phi_a.size(); i++) out.push_back(Phi_a[i]);
    for (int i = 0; i < Phi_b.size(); i++) out.push_back(Phi_b[i]);
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
    for (int i = 0; i < Phi.size(); i++) {
        if (Phi[i].spin() == spin) {
            out.push_back(Phi[i]);
        } else {
            tmp.push_back(Phi[i]);
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
 * @param n_orbs: number of orbitals to save
 *
 * The given file name (e.g. "phi") will be appended with orbital number ("phi_0").
 * Produces separate files for meta data ("phi_0.meta"), real ("phi_0_re.tree") and
 * imaginary ("phi_0_im.tree") parts. Negative n_orbs means that all orbitals in the
 * vector are saved.
 */
void orbital::save_orbitals(OrbitalVector &Phi, const std::string &file, int n_orbs) {
    if (n_orbs < 0) n_orbs = Phi.size();
    if (n_orbs > Phi.size()) MSG_ERROR("Index out of bounds");
    for (int i = 0; i < n_orbs; i++) {
        if (not mpi::my_orb(Phi[i])) continue; //only save own orbitals
        std::stringstream orbname;
        orbname << file << "_" << i;
        Phi[i].saveOrbital(orbname.str());
    }
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
    OrbitalVector Phi;
    for (int i = 0; true; i++) {
        if (n_orbs > 0 and i >= n_orbs) break;
        Orbital phi_i;
        std::stringstream orbname;
        orbname << file << "_" << i;
        phi_i.loadOrbital(orbname.str());
        phi_i.setRankId(mpi::orb_rank);
        if (phi_i.hasReal() or phi_i.hasImag()) {
            phi_i.setRankId(i % mpi::orb_size);
            Phi.push_back(phi_i);
            if (i % mpi::orb_size != mpi::orb_rank) phi_i.free();
        } else {
            break;
        }
    }
    //distribute errors
    DoubleVector errors = DoubleVector::Zero(Phi.size());
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) errors(i) = Phi[i].error();
    }
    mpi::allreduce_vector(errors, mpi::comm_orb);
    orbital::set_errors(Phi, errors);
    return Phi;
}

/** @brief Frees each orbital in the vector
 *
 * Leaves an empty vector. Orbitals are freed.
 *
 */
void orbital::free(OrbitalVector &Phi) {
    for (auto &phi_i : Phi) phi_i.free();
    Phi.clear();
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
void orbital::orthogonalize(Orbital &phi, Orbital psi) {
    ComplexDouble overlap = orbital::dot(psi, phi);
    double sq_norm = psi.squaredNorm();
    if (std::abs(overlap) > mrcpp::MachineZero) phi.add(-1.0 * overlap / sq_norm, psi);
}

/** @brief Gram-Schmidt orthogonalize orbitals within the set */
void orbital::orthogonalize(OrbitalVector &Phi) {
    mpi::free_foreign(Phi);
    for (int i = 0; i < Phi.size(); i++) {
        for (int j = 0; j < i; j++) {
            int tag = 7632 * i + j;
            int src = Phi[j].rankID();
            int dst = Phi[i].rankID();
            if (src == dst) {
                if (mpi::my_orb(Phi[i])) orbital::orthogonalize(Phi[i], Phi[j]);
            } else {
                if (mpi::my_orb(Phi[i])) {
                    mpi::recv_function(Phi[j], src, tag, mpi::comm_orb);
                    orbital::orthogonalize(Phi[i], Phi[j]);
                    Phi[j].free();
                }
                if (mpi::my_orb(Phi[j])) mpi::send_function(Phi[j], dst, tag, mpi::comm_orb);
            }
        }
    }
}

/** @brief Orthogonalize the Phi orbital against all orbitals in Psi */
void orbital::orthogonalize(OrbitalVector &Phi, OrbitalVector &Psi) {
    // Get all output orbitals belonging to this MPI
    OrbitalChunk myPhi = mpi::get_my_chunk(Phi);

    // Orthogonalize MY orbitals with ALL input orbitals
    OrbitalIterator iter(Psi, false);
    while (iter.next()) {
        for (int i = 0; i < iter.get_size(); i++) {
            Orbital &psi_i = iter.orbital(i);
            for (int j = 0; j < myPhi.size(); j++) {
                Orbital &phi_j = std::get<1>(myPhi[j]);
                orbital::orthogonalize(phi_j, psi_i);
            }
        }
    }
}

ComplexMatrix orbital::calc_overlap_matrix(OrbitalVector &BraKet) {
    ComplexMatrix S = ComplexMatrix::Zero(BraKet.size(), BraKet.size());

    // Get all ket orbitals belonging to this MPI
    OrbitalChunk myKet = mpi::get_my_chunk(BraKet);

    // Receive ALL orbitals on the bra side, use only MY orbitals on the ket side
    // Computes the FULL columns associated with MY orbitals on the ket side
    OrbitalIterator iter(BraKet, true); // use symmetry
    while (iter.next()) {
        for (int i = 0; i < iter.get_size(); i++) {
            int idx_i = iter.idx(i);
            Orbital &bra_i = iter.orbital(i);
            for (int j = 0; j < myKet.size(); j++) {
                int idx_j = std::get<0>(myKet[j]);
                Orbital &ket_j = std::get<1>(myKet[j]);
                if (mpi::my_orb(bra_i) and idx_j > idx_i) continue;
                if (mpi::my_unique_orb(ket_j) or mpi::orb_rank == 0) {
                    S(idx_i, idx_j) = orbital::dot(bra_i, ket_j);
                    S(idx_j, idx_i) = std::conj(S(idx_i, idx_j));
                }
            }
        }
    }
    // Assumes all MPIs have (only) computed their own columns of the matrix
    mpi::allreduce_matrix(S, mpi::comm_orb);
    return S;
}

/** @brief Compute the overlap matrix S_ij = <bra_i|ket_j>
 *
 * MPI: Each rank will compute the full columns related to their
 *      orbitals in the ket vector. The bra orbitals are communicated
 *      one rank at the time (all orbitals belonging to a given rank
 *      is communicated at the same time). This algorithm sets NO
 *      restrictions on the distributions of the bra or ket orbitals
 *      among the available ranks. After the columns have been computed,
 *      the full matrix is allreduced, e.i. all MPIs will have the full
 *      matrix at exit.
 *
 */
ComplexMatrix orbital::calc_overlap_matrix(OrbitalVector &Bra, OrbitalVector &Ket) {
    ComplexMatrix S = ComplexMatrix::Zero(Bra.size(), Ket.size());

    // Get all ket orbitals belonging to this MPI
    OrbitalChunk myKet = mpi::get_my_chunk(Ket);

    // Receive ALL orbitals on the bra side, use only MY orbitals on the ket side
    // Computes the FULL columns associated with MY orbitals on the ket side
    OrbitalIterator iter(Bra);
    while (iter.next()) {
        for (int i = 0; i < iter.get_size(); i++) {
            int idx_i = iter.idx(i);
            Orbital &bra_i = iter.orbital(i);
            for (int j = 0; j < myKet.size(); j++) {
                int idx_j = std::get<0>(myKet[j]);
                Orbital &ket_j = std::get<1>(myKet[j]);
                if (mpi::my_unique_orb(ket_j) or mpi::grand_master()) S(idx_i, idx_j) = orbital::dot(bra_i, ket_j);
            }
        }
    }
    // Assumes all MPIs have (only) computed their own columns of the matrix
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
    Timer timer;
    printout(1, "Calculating Löwdin orthonormalization matrix      ");

    ComplexMatrix S_tilde = orbital::calc_overlap_matrix(Phi);
    ComplexMatrix S_m12 = math_utils::hermitian_matrix_pow(S_tilde, -1.0 / 2.0);

    timer.stop();
    println(1, timer.getWallTime());
    return S_m12;
}

/** @brief Minimize the spatial extension of orbitals, by orbital rotation
 *
 * @param Phi: orbitals to localize
 *
 * Minimizes \f$  \sum_{i=1,N}\langle i| {\bf R^2}  | i \rangle - \langle i| {\bf R}| i \rangle^2 \f$
 * which is equivalent to maximizing \f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 *
 * The resulting transformation includes the orthonormalization of the orbitals.
 * Orbitals are rotated in place, and the transformation matrix is returned.
 */
ComplexMatrix orbital::localize(double prec, OrbitalVector &Phi) {
    Printer::printHeader(0, "Localizing orbitals");
    Timer timer;

    ComplexMatrix U;
    int n_it = 0;
    if (Phi.size() > 1) {
        Timer rmat;
        RRMaximizer rr(prec, Phi);
        rmat.stop();
        Printer::printDouble(0, "Computing position matrices", rmat.getWallTime(), 5);
        Timer rr_t;
        n_it = rr.maximize();
        rr_t.stop();
        Printer::printDouble(0, "Computing Foster-Boys matrix", rr_t.getWallTime(), 5);

        if (n_it > 0) {
            println(0, " Converged after iteration   " << std::setw(30) << n_it);
            U = rr.getTotalU().transpose().cast<ComplexDouble>();
        } else {
            println(0, " Foster-Boys localization did not converge!");
        }
    } else {
        println(0, " Cannot localize less than two orbitals");
    }

    if (n_it <= 0) {
        Timer orth_t;
        U = orbital::calc_lowdin_matrix(Phi);
        orth_t.stop();
        Printer::printDouble(0, "Computing Lowdin matrix", orth_t.getWallTime(), 5);
    }

    Timer rot_t;
    OrbitalVector Psi = orbital::rotate(U, Phi, prec);
    orbital::free(Phi);
    Phi = Psi;
    rot_t.stop();
    Printer::printDouble(0, "Rotating orbitals", rot_t.getWallTime(), 5);

    timer.stop();
    Printer::printFooter(0, timer, 2);
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
    Printer::printHeader(0, "Digonalizing Fock matrix");
    Timer timer;

    Timer orth_t;
    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);
    F = S_m12.transpose() * F * S_m12;
    orth_t.stop();
    Printer::printDouble(0, "Computing Lowdin matrix", orth_t.getWallTime(), 5);

    Timer diag_t;
    ComplexMatrix U = ComplexMatrix::Zero(F.rows(), F.cols());
    int np = orbital::size_paired(Phi);
    int na = orbital::size_alpha(Phi);
    int nb = orbital::size_beta(Phi);
    if (np > 0) math_utils::diagonalize_block(F, U, 0, np);
    if (na > 0) math_utils::diagonalize_block(F, U, np, na);
    if (nb > 0) math_utils::diagonalize_block(F, U, np + na, nb);
    U = U * S_m12;
    diag_t.stop();
    Printer::printDouble(0, "Diagonalizing matrix", diag_t.getWallTime(), 5);

    Timer rot_t;
    OrbitalVector Psi = orbital::rotate(U, Phi, prec);
    orbital::free(Phi);
    Phi = Psi;
    rot_t.stop();
    Printer::printDouble(0, "Rotating orbitals", rot_t.getWallTime(), 5);

    timer.stop();
    Printer::printFooter(0, timer, 2);
    return U;
}

/** @brief Perform the Löwdin orthonormalization
 *
 * @param Phi: orbitals to orthonormalize
 *
 * Orthonormalizes the orbitals by multiplication of the Löwdin matrix S^(-1/2).
 * Orbitals are rotated in place, and the transformation matrix is returned.
 */
ComplexMatrix orbital::orthonormalize(double prec, OrbitalVector &Phi) {
    ComplexMatrix U = orbital::calc_lowdin_matrix(Phi);
    OrbitalVector Psi = orbital::rotate(U, Phi, prec);
    orbital::free(Phi);
    Phi = Psi;
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
        if (phi_i.occ() == 1) nDoubly++;
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

/** @brief Returns a vector containing the orbital errors */
DoubleVector orbital::get_errors(const OrbitalVector &Phi) {
    int nOrbs = Phi.size();
    DoubleVector errors = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) errors(i) = Phi[i].error();
    return errors;
}

/** @brief Assign errors to each orbital.
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void orbital::set_errors(OrbitalVector &Phi, const DoubleVector &errors) {
    if (Phi.size() != errors.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < Phi.size(); i++) Phi[i].setError(errors(i));
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

/** @brief Returns a vector containing the orbital occupancies */
IntVector orbital::get_occupancies(const OrbitalVector &Phi) {
    int nOrbs = Phi.size();
    IntVector occ = IntVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) occ(i) = Phi[i].occ();
    return occ;
}

/** @brief Assigns spin to each orbital
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void orbital::set_occupancies(OrbitalVector &Phi, const IntVector &occ) {
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
    return norms;
}

/** @brief Returns a vector containing the orbital norms */
DoubleVector orbital::get_norms(const OrbitalVector &Phi) {
    int nOrbs = Phi.size();
    DoubleVector norms = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        if (mpi::my_orb(Phi[i])) norms(i) = Phi[i].norm();
    }
    return norms;
}

void orbital::print(const OrbitalVector &Phi) {
    Printer::setScientific();
    printout(0, "============================================================\n");
    printout(0, " OrbitalVector:");
    printout(0, std::setw(4) << Phi.size() << " orbitals  ");
    printout(0, std::setw(4) << size_occupied(Phi) << " occupied  ");
    printout(0, std::setw(4) << get_electron_number(Phi) << " electrons\n");
    printout(0, "------------------------------------------------------------\n");
    printout(0, "   n  RankID           Norm          Spin Occ      Error    \n");
    printout(0, "------------------------------------------------------------\n");
    for (int i = 0; i < Phi.size(); i++) println(0, std::setw(4) << i << Phi[i]);
    printout(0, "============================================================\n\n\n");
}

} //namespace mrchem
