/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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
#include "HelmholtzVector.h"
#include "OrbitalOptimizer.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/two_electron/FockOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief constructor
 *
 * @param h: Helmholtz operators
 * @param k: Iterative accelerator
 *
 * SCF solver will NOT take ownership of the HelmholtzVector or the Accelerator,
 * so the original objects must be taken care of externally (do not delete until
 * SCF goes out of scope). Fock matrix, FockOperator and OrbitalVector are not
 * initialized at this stage, so the SCF solver needs to be "setup()" before
 * "optimize()".
 */
OrbitalOptimizer::OrbitalOptimizer(Accelerator *k)
        : GroundStateSolver()
        , kain(k) {}

/** @brief Prepare solver for optimization
 *
 * @param F: FockOperator defining the SCF equations
 * @param Phi: Orbitals to optimize
 * @param F_mat: Fock matrix
 *
 * SCF solver will NOT take ownership of the input, so these objects must be taken
 * care of externally (do not delete until SCF goes out of scope).
 */
void OrbitalOptimizer::setup(FockOperator &F, OrbitalVector &Phi, ComplexMatrix &F_mat) {
    this->fMat_n = &F_mat;
    this->fOper_n = &F;
    this->orbitals_n = &Phi;
}

/** @brief Clear solver after optimization
 *
 * Clear pointers that was set during setup, and reset the precision parameter
 * (only the current precision orbPrec[0], not the boundary values orbPrec[1,2]).
 * Solver can be re-used after another setup.
 */
void OrbitalOptimizer::clear() {
    this->fMat_n = nullptr;
    this->fOper_n = nullptr;
    this->orbitals_n = nullptr;
    if (useKAIN()) this->kain->clear();
    resetPrecision();
}

/** @brief Run orbital optimization
 *
 * Optimize orbitals until convergence thresholds are met. This algorithm computes
 * the Fock matrix explicitly using the kinetic energy operator, and uses a KAIN
 * accelerator to improve convergence. Diagonalization or localization may be performed
 * during the SCF iterations. Main points of the algorithm:
 *
 * Pre SCF: setup Fock operator and compute Fock matrix
 *
 *  1) Diagonalize/localize orbitals
 *  2) Compute current SCF energy
 *  3) Apply Helmholtz operator on all orbitals
 *  4) Orthonormalize orbitals (Löwdin)
 *  5) Compute orbital updates
 *  6) Compute KAIN update
 *  7) Compute errors and check for convergence
 *  8) Add orbital updates
 *  9) Orthonormalize orbitals (Löwdin)
 * 10) Setup Fock operator
 * 11) Compute Fock matrix
 *
 * Post SCF: diagonalize/localize orbitals
 */
bool OrbitalOptimizer::optimize() {
    ComplexMatrix &F_mat = *this->fMat_n;
    OrbitalVector &Phi_n = *this->orbitals_n;
    FockOperator &F = *this->fOper_n;
    RankZeroTensorOperator &V = F.potential();

    double orb_prec = this->orbPrec[0];
    double err_o = orbital::get_errors(Phi_n).maxCoeff();
    double err_t = 1.0;
    double err_p = 1.0;

    F.setup(orb_prec);
    F_mat = F(Phi_n, Phi_n);

    int nIter = 0;
    bool converged = false;
    while (nIter++ < this->maxIter or this->maxIter < 0) {
        // Initialize SCF cycle
        Timer timer;
        printCycleHeader(nIter);
        orb_prec = adjustPrecision(err_o);

        // Rotate orbitals
        if (needLocalization(nIter)) {
            ComplexMatrix U_mat = orbital::localize(orb_prec, Phi_n);
            F.rotate(U_mat);
            F_mat = U_mat * F_mat * U_mat.adjoint();
            if (useKAIN()) this->kain->clear();
        } else if (needDiagonalization(nIter)) {
            ComplexMatrix U_mat = orbital::diagonalize(orb_prec, Phi_n, F_mat);
            F.rotate(U_mat);
            if (useKAIN()) this->kain->clear();
        }
        // Compute electronic energy
        double E = calcProperty();
        this->property.push_back(E);

        // Apply Helmholtz operator
        HelmholtzVector H(orb_prec, F_mat.real().diagonal());
        OrbitalVector Psi = H.rotate(F_mat, Phi_n);
        OrbitalVector Phi_np1 = H(V, Phi_n, Psi);
        Psi.clear();
        F.clear();

        ComplexMatrix U_mat = orbital::orthonormalize(orb_prec, Phi_np1);
        F_mat = U_mat * F_mat * U_mat.adjoint();

        // Compute orbital updates
        OrbitalVector dPhi_n = orbital::add(1.0, Phi_np1, -1.0, Phi_n);
        Phi_np1.clear();

        // Employ KAIN accelerator
        if (useKAIN()) this->kain->accelerate(orb_prec, Phi_n, dPhi_n);

        // Compute errors
        DoubleVector errors = orbital::get_norms(dPhi_n);

        err_o = errors.maxCoeff();
        err_t = errors.norm();
        err_p = calcPropertyError();
        this->orbError.push_back(err_t);
        converged = checkConvergence(err_o, err_p);

        // Update orbitals
        Phi_n = orbital::add(1.0, Phi_n, 1.0, dPhi_n);
        dPhi_n.clear();

        orbital::orthonormalize(orb_prec, Phi_n);
        orbital::set_errors(Phi_n, errors);

        // Compute Fock matrix
        F.setup(orb_prec);
        F_mat = F(Phi_n, Phi_n);

        // Finalize SCF cycle
        timer.stop();
        printOrbitals(F_mat.real().diagonal(), Phi_n, 0);
        printProperty();
        printCycleFooter(timer.getWallTime());

        if (converged) break;
    }

    if (useKAIN()) this->kain->clear();
    F.clear();

    printConvergence(converged);
    return converged;
}

} // namespace mrchem
