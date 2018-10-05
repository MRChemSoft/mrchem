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

#include "Accelerator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

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
    while (this->orbitals.size() > 0) {
        orbital::free(this->orbitals[0]);
        this->orbitals.pop_front();
    }
    while (this->dOrbitals.size() > 0) {
        orbital::free(this->dOrbitals[0]);
        this->dOrbitals.pop_front();
    }
    while (this->fock.size() > 0) { this->fock.pop_front(); }
    while (this->dFock.size() > 0) { this->dFock.pop_front(); }
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
    Timer timer;
    int nOrbs = this->orbitals.size() - 1;
    int nFock = this->fock.size() - 1;
    if (all) {
        nOrbs += 1;
        nFock += 1;
    }
    if (nOrbs <= 0) { return; }
    for (int i = 0; i < nOrbs; i++) {
        OrbitalVector &Phi = this->orbitals[i];
        OrbitalVector tmp = orbital::rotate(U, Phi);
        orbital::free(Phi);
        Phi = tmp;

        OrbitalVector &dPhi = this->dOrbitals[i];
        OrbitalVector dTmp = orbital::rotate(U, dPhi);
        orbital::free(dPhi);
        dPhi = dTmp;
    }
    for (int i = 0; i < nFock; i++) {
        ComplexMatrix &F = this->fock[i];
        ComplexMatrix &dF = this->dFock[i];
        F = U * F * U.transpose();
        dF = U * dF * U.transpose();
    }

    timer.stop();
    double t = timer.getWallTime();
    Printer::printDouble(0, "Rotating iterative subspace", t, 5);
}

/** @brief Update iterative history with the latest orbitals and updates
 *
 * @param Phi: Next set ov orbitals
 * @param dPhi: Next set of orbital updates
 * @param F: Next Fock matrix
 * @param dF: Next Fock matrix update
 *
 * The new orbitals are moved from the input sets and into the newly
 * constructed local sets. This means that the input sets contain a
 * bunch of zero pointers after this routine. To get the orbitals
 * back, use calcUpdates or getOrbitals. If F and dF are given as
 * input, the matrices are included in the subspace. If the length
 * of the history exceed maxHistory the oldest orbitals are discarded.
 */
// clang-format off
void Accelerator::push_back(OrbitalVector &Phi,
                            OrbitalVector &dPhi,
                            ComplexMatrix *F,
                            ComplexMatrix *dF) {
    // clang-format on
    Timer timer;
    int nHistory = this->orbitals.size();
    bool historyIsFull = false;
    if (nHistory >= this->maxHistory) { historyIsFull = true; }
    if (F != 0) {
        if (dF == 0) { MSG_ERROR("Need to give both F and dF"); }
        if (this->fock.size() != nHistory) { MSG_ERROR("Size mismatch orbitals vs matrices"); }
    }

    if (this->orbitals.size() > this->maxHistory - 1) {
        orbital::free(this->orbitals[0]);
        this->orbitals.pop_front();
    }
    if (this->dOrbitals.size() > this->maxHistory - 1) {
        orbital::free(this->dOrbitals[0]);
        this->dOrbitals.pop_front();
    }
    if (this->fock.size() > this->maxHistory - 1) { this->fock.pop_front(); }
    if (this->dFock.size() > this->maxHistory - 1) { this->dFock.pop_front(); }
    if (not verifyOverlap(Phi)) {
        println(0, " Clearing accelerator");
        this->clear();
    }

    this->orbitals.push_back(Phi);
    this->dOrbitals.push_back(dPhi);
    if (F != 0) this->fock.push_back(*F);
    if (dF != 0) this->dFock.push_back(*dF);

    Phi.clear();
    dPhi.clear();

    timer.stop();
    double t = timer.getWallTime();
    Printer::printDouble(0, "Push back orbitals", t, 5);
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
    IntVector out = IntVector::Zero(nOrbs);
    int nHistory = this->orbitals.size() - 1;
    if (nHistory > 0) {
        for (int i = 0; i < nOrbs; i++) {
            Orbital &phi_i = Phi[i];
            if (mpi::my_orb(phi_i)) {
                Orbital &last_i = this->orbitals[nHistory][i];
                if (not mpi::my_orb(last_i)) MSG_FATAL("MPI rank mismatch");
                double sqNorm = phi_i.squaredNorm();
                ComplexDouble overlap = orbital::dot(phi_i, last_i);
                if (std::abs(overlap) < 0.5 * sqNorm) {
                    Printer::printDouble(0, "Overlap not verified ", std::abs(overlap));
                    out(i) = 1;
                }
            }
        }
    }
    mpi::allreduce_vector(out, mpi::comm_orb);

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
 * Solves the linear system of equations \f$ Ac = b \f$ that defines the
 * coefficients \f$ c \f$ of the next step
 * copies the resulting orbitals into the output sets. If length of history
 * is less than minHistory, the latest orbitals are copied directly into
 * the output sets without solving the linear problem. Existing input
 * orbitals and updates are replaced by new ones.
 */
void Accelerator::accelerate(double prec,
                             OrbitalVector &Phi,
                             OrbitalVector &dPhi,
                             ComplexMatrix *F,
                             ComplexMatrix *dF) {
    Printer::printHeader(0, "Iterative subspace accelerator");

    Timer timer;
    this->push_back(Phi, dPhi, F, dF);

    int nHistory = this->orbitals.size() - 1;
    if (nHistory <= this->minHistory) {
        copyOrbitals(Phi);
        copyOrbitalUpdates(dPhi);
    } else {
        setupLinearSystem();
        solveLinearSystem();
        expandSolution(prec, Phi, dPhi, F, dF);
        clearLinearSystem();
    }
    timer.stop();
    Printer::printFooter(0, timer, 2);
}

/** @brief Solve the linear system of equations
 *
 * Solves the linear system of equations \f$ Ac = b \f$ that defines the
 * coefficients \f$ c \f$ of the next step.
 */
void Accelerator::solveLinearSystem() {
    Timer timer;
    this->c.clear();
    int N = this->A.size();
    for (int n = 0; n < N; n++) {
        DoubleVector tmpC = this->A[n].colPivHouseholderQr().solve(this->b[n]);
        this->c.push_back(tmpC);
    }
    timer.stop();
    double t = timer.getWallTime();
    Printer::printDouble(0, "Solve linear system", t, 5);
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
    Timer timer;
    int totHistory = this->orbitals.size();
    if (nHistory >= totHistory or nHistory < 0) { MSG_FATAL("Requested orbitals unavailable"); }
    int n = totHistory - 1 - nHistory;
    Phi = orbital::deep_copy(this->orbitals[n]);

    timer.stop();
    double t = timer.getWallTime();
    Printer::printDouble(0, "Copy orbitals", t, 5);
}

/** @brief Copy orbital updates at the given point in history into the given set.
 *
 * @param dPhi: OrbitalVector to copy into
 * @param nHistory: point in history to copy
 *
 * Input set is overwritten if it contains orbitals. Counts backwards,
 * zero history input returns latest orbital set. */
void Accelerator::copyOrbitalUpdates(OrbitalVector &dPhi, int nHistory) {
    Timer timer;
    int totHistory = this->dOrbitals.size();
    if (nHistory >= totHistory or nHistory < 0) { MSG_FATAL("Requested orbitals unavailable"); }
    int n = totHistory - 1 - nHistory;
    dPhi = orbital::deep_copy(this->dOrbitals[n]);

    timer.stop();
    double t = timer.getWallTime();
    Printer::printDouble(0, "Copy orbital updates", t, 5);
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
    //        MSG_FATAL("Requested orbitals unavailable");
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
    //        MSG_FATAL("Requested orbitals unavailable");
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
// clang-format off
void Accelerator::sortLinearSystem(std::vector<DoubleMatrix> &A_matrices,
                                   std::vector<DoubleVector> &b_vectors) {
    // clang-format on
    int nOrbs = b_vectors.size();
    int nHist = b_vectors[0].size();

    if (this->sepOrbitals) {
        for (int i = 0; i < nOrbs; i++) {
            DoubleMatrix tmpA(A_matrices[i]);
            DoubleVector tmpB(b_vectors[i]);
            this->A.push_back(tmpA);
            this->b.push_back(tmpB);
        }
    } else {
        DoubleMatrix tmpA(nHist, nHist);
        DoubleVector tmpB(nHist);
        tmpA.setZero();
        tmpB.setZero();
        for (int i = 0; i < nOrbs; i++) {
            tmpA += A_matrices[i];
            tmpB += b_vectors[i];
        }
        this->A.push_back(tmpA);
        this->b.push_back(tmpB);
    }
}

/** @brief Prints the number of trees and nodes kept in the iterative history */
int Accelerator::printTreeSizes() const {
    int totNodes = 0;
    int totTrees = 0;
    int nHistory = this->orbitals.size();
    for (int i = 0; i < nHistory; i++) {
        int nOrbitals = this->orbitals[i].size();
        for (int j = 0; j < nOrbitals; j++) { totNodes += this->orbitals[i][j].getNNodes(); }
        totTrees += nOrbitals;
    }
    for (int i = 0; i < nHistory; i++) {
        int nOrbitals = this->dOrbitals[i].size();
        for (int j = 0; j < nOrbitals; j++) { totNodes += this->dOrbitals[i][j].getNNodes(); }
        totTrees += nOrbitals;
    }
    println(0, " Accelerator       " << std::setw(15) << totTrees << std::setw(25) << totNodes);
    return totNodes;
}

} // namespace mrchem
