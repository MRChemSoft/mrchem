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

#include "Accelerator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief constructor
 *
 * @param max: max length of history
 * @param min: min length of history
 * @param sep: solve separate eqations for each orbital
 */
Accelerator::Accelerator(int max, int min, bool sep)
        : minHistory(min)
        , maxHistory(max)
        , sepOrbitals(sep) {
    if (this->minHistory < 0) this->minHistory = 0;
    if (this->maxHistory < this->minHistory) MSG_ERROR("Invalid argument");
}

/** @brief Delete all orbitals and Fock matrices in history.
 *
 * Leaves the accelerator in the same state as after construction,
 * ready for use in the next optimization.
 */
void Accelerator::clear() {
    while (this->orbitals.size() > 0) this->orbitals.pop_front();
    while (this->dOrbitals.size() > 0) this->dOrbitals.pop_front();
    while (this->fock.size() > 0) this->fock.pop_front();
    while (this->dFock.size() > 0) this->dFock.pop_front();
    clearLinearSystem();
}

/** @brief Delete the matrices and vectors used to compute the next step.
 *
 * The accelerator is now ready to move to the next iteration.
 */
void Accelerator::clearLinearSystem() {
    this->A.clear();
    this->b.clear();
    this->c.clear();
}

/** @brief Rotate iterative history.
 *
 * @param U: rotation matrix
 * @param all: rotate ALL orbitals in the history
 *
 * To keep phases and orbital ordering consistent one can apply
 * the latest orbital rotation to the entire orbital history.
 * Not recommended as rotations are expensive. Instead one should
 * clear history and start over. Option to rotate the last orbital
 * set or not.
 */
void Accelerator::rotate(const ComplexMatrix &U, bool all) {
    Timer t_tot;
    int nOrbs = this->orbitals.size() - 1;
    int nFock = this->fock.size() - 1;
    if (all) {
        nOrbs += 1;
        nFock += 1;
    }
    if (nOrbs <= 0) { return; }
    for (int i = 0; i < nOrbs; i++) {
        auto &Phi = this->orbitals[i];
        mrcpp::mpifuncvec::rotate(Phi, U);

        auto &dPhi = this->dOrbitals[i];
        mrcpp::mpifuncvec::rotate(dPhi, U);
    }
    for (int i = 0; i < nFock; i++) {
        auto &F = this->fock[i];
        auto &dF = this->dFock[i];
        F = U.adjoint() * F * U;
        dF = U.adjoint() * dF * U;
    }
    mrcpp::print::time(this->pl + 2, "Rotating iterative subspace", t_tot);
}

/** @brief Update iterative history with the latest orbitals and updates
 *
 * @param Phi: Next set of orbitals
 * @param dPhi: Next set of orbital updates
 * @param F: Next Fock matrix
 * @param dF: Next Fock matrix update
 *
 * The new orbitals are deep copied from the input sets and into the
 * accelerator history. The input sets are left unchanged. If F and
 * dF are given as input, the matrices are included in the subspace.
 * If the length of the history exceed maxHistory the oldest orbitals
 * are discarded.
 */
void Accelerator::push_back(OrbitalVector &Phi, OrbitalVector &dPhi, ComplexMatrix *F, ComplexMatrix *dF) {
    Timer t_tot;
    int nHistory = this->orbitals.size();
    if (F != nullptr) {
        if (dF == nullptr) MSG_ERROR("Need to give both F and dF");
        if (this->fock.size() != nHistory) MSG_ERROR("Size mismatch orbitals vs matrices");
    }
    auto historyIsFull = (nHistory >= this->maxHistory);
    if (historyIsFull and this->orbitals.size() > 0) this->orbitals.pop_front();
    if (historyIsFull and this->dOrbitals.size() > 0) this->dOrbitals.pop_front();
    if (historyIsFull and this->fock.size() > 0) this->fock.pop_front();
    if (historyIsFull and this->dFock.size() > 0) this->dFock.pop_front();

    if (not verifyOverlap(Phi)) {
        println(this->pl + 2, " Clearing accelerator");
        this->clear();
    }

    this->orbitals.push_back(orbital::deep_copy(Phi));
    this->dOrbitals.push_back(orbital::deep_copy(dPhi));
    if (F != nullptr) this->fock.push_back(*F);
    if (dF != nullptr) this->dFock.push_back(*dF);

    mrcpp::print::time(this->pl + 2, "Push back orbitals", t_tot);
}

/** @brief Verify that the orbital overlap between the two last iterations is positive.
 *
 * @param Phi: Next set of orbitals
 *
 * If the overlap is negative a sign change of the corresponding
 * orbitals in the entire history is required (unless the history
 * is cleared).
 */
bool Accelerator::verifyOverlap(OrbitalVector &Phi) {
    int nOrbs = Phi.size();
    int nHistory = this->orbitals.size() - 1;
    auto out = IntVector::Zero(nOrbs).eval();
    if (nHistory > 0) {
        for (int i = 0; i < nOrbs; i++) {
            auto &phi_i = Phi[i];
            if (mrcpp::mpi::my_orb(phi_i)) {
                auto &last_i = this->orbitals[nHistory][i];
                if (not mrcpp::mpi::my_orb(last_i)) MSG_ABORT("MPI rank mismatch");
                auto sqNorm = phi_i.squaredNorm();
                auto overlap = orbital::dot(phi_i, last_i);
                if (std::abs(overlap) < 0.5 * sqNorm) {
                    mrcpp::print::value(this->pl + 2, "Overlap not verified ", std::abs(overlap));
                    out(i) = 1;
                }
            }
        }
    }
    mrcpp::mpi::allreduce_vector(out, mrcpp::mpi::comm_wrk);

    return (out.sum() < 1) ? true : false;
}

/** @brief Calculates the new orbitals and updates based on history information
 *
 * @param prec: Precision used in arithmetic operations
 * @param Phi: Orbitals to accelerate (in/out)
 * @param dPhi: Orbital updates to accelerate (in/out)
 * @param F: Fock matrix to accelerate (in/out)
 * @param dF: Fock matrix update to accelerate (in/out)
 *
 * This routine will make a local copy of the input orbitals/matrices
 * and keep the history. It then solves a linear system of equations
 * \f$ Ac = b \f$ that defines the coefficients \f$ c \f$ of the next
 * step. If length of history is _smaller_ than minHistory, the input
 * orbitals/matrices remain _unchanged_, otherwise they are _overwritten_
 * by new ones. If the length of history is _larger_ than maxHistory, the
 * oldest iteration is discarded.
 */
void Accelerator::accelerate(double prec, OrbitalVector &Phi, OrbitalVector &dPhi, ComplexMatrix *F, ComplexMatrix *dF) {
    if (this->maxHistory < 1) return;

    Timer t_tot;
    auto plevel = Printer::getPrintLevel(); // global print level
    mrcpp::print::header(this->pl + 2, "Iterative subspace accelerator");

    // Deep copy into history
    this->push_back(Phi, dPhi, F, dF);

    int nHistory = this->orbitals.size() - 1;
    if (nHistory > this->minHistory) {
        setupLinearSystem();
        solveLinearSystem();
        // Overwrites (Phi, dPhi, F, dF) with new guess
        expandSolution(prec, Phi, dPhi, F, dF);
        clearLinearSystem();
    }
    printSizeNodes();
    mrcpp::print::footer(this->pl + 2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(this->pl + 1, "Iterative subspace accelerator", t_tot);
}

/** @brief Solve the linear system of equations
 *
 * Solves the linear system of equations \f$ Ac = b \f$ that defines the
 * coefficients \f$ c \f$ of the next step.
 */
void Accelerator::solveLinearSystem() {
    Timer t_tot;
    this->c.clear();
    for (int n = 0; n < this->A.size(); n++) {
        auto tmpC = this->A[n].colPivHouseholderQr().solve(this->b[n]).eval();
        this->c.push_back(tmpC);
    }
    mrcpp::print::time(this->pl + 2, "Solve linear system", t_tot);
}

/** @brief Copy orbitals at the given point in history into the given set.
 *
 * @param Phi: OrbitalVector to copy into
 * @param nHistory: point in history to copy
 *
 * Input set is overwritten if it contains orbitals. Counts backwards,
 * zero history input returns latest orbital set.
 */
void Accelerator::copyOrbitals(OrbitalVector &Phi, int nHistory) {
    Timer t_tot;
    int totHistory = this->orbitals.size();
    if (nHistory >= totHistory or nHistory < 0) MSG_ABORT("Requested orbitals unavailable");
    int n = totHistory - 1 - nHistory;
    Phi = orbital::deep_copy(this->orbitals[n]);
    mrcpp::print::time(this->pl + 2, "Copy orbitals", t_tot);
}

/** @brief Copy orbital updates at the given point in history into the given set.
 *
 * @param dPhi: OrbitalVector to copy into
 * @param nHistory: point in history to copy
 *
 * Input set is overwritten if it contains orbitals. Counts backwards,
 * zero history input returns latest orbital set. */
void Accelerator::copyOrbitalUpdates(OrbitalVector &dPhi, int nHistory) {
    Timer t_tot;
    int totHistory = this->dOrbitals.size();
    if (nHistory >= totHistory or nHistory < 0) MSG_ABORT("Requested orbitals unavailable");
    int n = totHistory - 1 - nHistory;
    dPhi = orbital::deep_copy(this->dOrbitals[n]);
    mrcpp::print::time(this->pl + 2, "Copy orbital updates", t_tot);
}

/** @brief Replaces the orbital set from a given point in history.
 *
 * @param Phi: OrbitalVector to copy into
 * @param nHistory: point in history to copy
 *
 * Deletes the old orbital set and copies the new. Counts backwards,
 * zero input returns latest orbital set. */
void Accelerator::replaceOrbitals(OrbitalVector &Phi, int nHistory) {
    NOT_IMPLEMENTED_ABORT;
    //    int totHistory = this->orbitals.size();
    //    if (nHistory >= totHistory or nHistory < 0) {
    //        MSG_ABORT("Requested orbitals unavailable");
    //    }
    //    int n = totHistory - 1 - nHistory;

    //    string oldName = this->orbitals[n]->getName();
    //    OrbitalSet *orbSet = new OrbitalSet(oldName, orbs);
    //    *orbSet = orbs;
    //    delete this->orbitals[n];
    //    this->orbitals[n] = orbSet;
}

/** @brief Replaces the orbital update set from a given point in history
 *
 * @param dPhi: OrbitalVector to copy into
 * @param nHistory: point in history to copy
 *
 * Deletes the old orbital set and copies the new. Counts backwards,
 * zero input returns latest orbital set. */
void Accelerator::replaceOrbitalUpdates(OrbitalVector &dPhi, int nHistory) {
    NOT_IMPLEMENTED_ABORT;
    //    int totHistory = this->dOrbitals.size();
    //    if (nHistory >= totHistory or nHistory < 0) {
    //        MSG_ABORT("Requested orbitals unavailable");
    //    }
    //    int n = totHistory - 1 - nHistory;

    //    string oldName = this->dOrbitals[n]->getName();
    //    OrbitalSet *orbSet = new OrbitalSet(oldName, dOrbs);
    //    *orbSet = dOrbs;
    //    delete this->dOrbitals[n];
    //    this->dOrbitals[n] = orbSet;
}

/** @brief Creates the final A matrices and b vectors.
 *
 * @param A_matrices: Matrices to collect
 * @param b_vectors: Vectors to collect
 *
 * Input vectors may include nOrbs + 1 entries, one extra from the Fock matrix.
 * If orbitals are not separated these are added up to one final A matrix and
 * b vector, otherwise all individual entries are kept.
 */
void Accelerator::sortLinearSystem(std::vector<ComplexMatrix> &A_matrices, std::vector<ComplexVector> &b_vectors) {
    if (this->sepOrbitals) {
        for (int i = 0; i < b_vectors.size(); i++) {
            auto tmpA(A_matrices[i]);
            auto tmpB(b_vectors[i]);
            this->A.push_back(tmpA);
            this->b.push_back(tmpB);
        }
    } else {
        auto tmpA(A_matrices[0]);
        auto tmpB(b_vectors[0]);
        tmpA.setZero();
        tmpB.setZero();
        for (int i = 0; i < b_vectors.size(); i++) {
            tmpA += A_matrices[i];
            tmpB += b_vectors[i];
        }
        this->A.push_back(tmpA);
        this->b.push_back(tmpB);
    }
}

/** @brief Prints the number of trees and nodes kept in the iterative history */
void Accelerator::printSizeNodes() const {
    int n = 0, m = 0;
    for (const auto &orbs_i : this->orbitals) {
        n += orbital::get_n_nodes(orbs_i);
        m += orbital::get_size_nodes(orbs_i);
    }
    int dn = 0, dm = 0;
    for (const auto &orbs_i : this->dOrbitals) {
        dn += orbital::get_n_nodes(orbs_i);
        dm += orbital::get_size_nodes(orbs_i);
    }
    mrcpp::print::separator(this->pl + 2, '-');
    mrcpp::print::tree(this->pl + 2, "Orbital sizes", n, m, 0.0);
    mrcpp::print::tree(this->pl + 2, "Update sizes", dn, dm, 0.0);
}

} // namespace mrchem
