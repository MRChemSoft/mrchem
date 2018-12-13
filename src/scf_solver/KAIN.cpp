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

#include "KAIN.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Calculates the A matrices and b vectors for each orbital individually.
 *
 * \f$ A_{ij} = \langle x^n - x^i | f(x^n) - f(x^j) \rangle \f$
 * \f$ b_{i} = \langle x^n - x^i | f(x^n) \rangle \f$
 *
 * If Fock matrix is included this is treated as an additional ''orbital''
 * and the return vectors have size nOrbs + 1. Frobenius inner product
 * used for the Fock matrix. If separateOrbitals is false the A's and b's
 * are later collected to single entities.
 */
void KAIN::setupLinearSystem() {
    Timer timer;
    int nHistory = this->orbitals.size() - 1;
    if (nHistory < 1) MSG_FATAL("Not enough history to setup system of equations");

    std::vector<DoubleMatrix> A_matrices;
    std::vector<DoubleVector> b_vectors;

    int nOrbitals = this->orbitals[nHistory].size();
    for (int n = 0; n < nOrbitals; n++) {
        DoubleMatrix orbA = DoubleMatrix::Zero(nHistory, nHistory);
        DoubleVector orbB = DoubleVector::Zero(nHistory);

        Orbital &phi_m = this->orbitals[nHistory][n];
        Orbital &fPhi_m = this->dOrbitals[nHistory][n];
        if (mpi::my_orb(phi_m)) {
            if (not mpi::my_orb(fPhi_m)) MSG_FATAL("MPI rank mismatch: fPhi_m");

            for (int i = 0; i < nHistory; i++) {
                Orbital &phi_i = this->orbitals[i][n];
                if (not mpi::my_orb(phi_i)) MSG_FATAL("MPI rank mismatch: phi_i");
                Orbital dPhi_im = phi_m.paramCopy();
                qmfunction::add(dPhi_im, 1.0, phi_i, -1.0, phi_m, -1.0);

                for (int j = 0; j < nHistory; j++) {
                    Orbital &fPhi_j = this->dOrbitals[j][n];
                    if (not mpi::my_orb(fPhi_j)) MSG_FATAL("MPI rank mismatch: fPhi_j");
                    Orbital dfPhi_jm = fPhi_m.paramCopy();
                    qmfunction::add(dfPhi_jm, 1.0, fPhi_j, -1.0, fPhi_m, -1.0);

                    // Ref. Harrisons KAIN paper the following has the wrong sign,
                    // but we define the updates (lowercase f) with opposite sign.
                    orbA(i, j) -= orbital::dot(dPhi_im, dfPhi_jm).real();
                }
                orbB(i) += orbital::dot(dPhi_im, fPhi_m).real();
            }
        }
        A_matrices.push_back(orbA);
        b_vectors.push_back(orbB);
    }

    for (int i = 0; i < nOrbitals; i++) {
        mpi::allreduce_matrix(A_matrices[i], mpi::comm_orb);
        mpi::allreduce_vector(b_vectors[i], mpi::comm_orb);
    }

    // Fock matrix is treated as a whole using the Frobenius inner product
    if (this->orbitals.size() == this->fock.size()) {
        const ComplexMatrix &X_m = this->fock[nHistory];
        const ComplexMatrix &fX_m = this->dFock[nHistory];

        ComplexMatrix fockA = ComplexMatrix::Zero(nHistory, nHistory);
        ComplexVector fockB = ComplexVector::Zero(nHistory);

        for (int i = 0; i < nHistory; i++) {
            const ComplexMatrix &X_i = this->fock[i];
            ComplexMatrix dX_im = X_i - X_m;
            for (int j = 0; j < nHistory; j++) {
                const ComplexMatrix &fX_j = this->dFock[j];
                ComplexMatrix dfX_jm = fX_j - fX_m;
                ComplexMatrix prod = dX_im.transpose() * dfX_jm;
                fockA(i, j) -= prod.trace();
            }
            ComplexMatrix prod = dX_im.transpose() * fX_m;
            fockB(i) += prod.trace();
        }
        A_matrices.push_back(fockA.real());
        b_vectors.push_back(fockB.real());
    }

    sortLinearSystem(A_matrices, b_vectors);

    timer.stop();
    double t = timer.getWallTime();
    Printer::printDouble(0, "Setup linear system", t, 5);
}

/** @brief Compute the next step for orbitals and orbital updates
 *
 * The next step \f$ \delta x^n \f$ is constructed from the solution
 * solution \f$ c \f$ of the linear problem \f$ Ac = b \f$ as:
 *
 * \f$ \delta x^n = f(x^n) + \sum_{j=1}^m c_j[(x^j-x^n)+(f(x^j)-f(x^n))]\f$
 */
void KAIN::expandSolution(double prec, OrbitalVector &Phi, OrbitalVector &dPhi, ComplexMatrix *F, ComplexMatrix *dF) {
    Timer timer;
    int nHistory = this->orbitals.size() - 1;
    int nOrbitals = this->orbitals[nHistory].size();

    // Orbitals are unchanged, updates will change
    Phi = orbital::deep_copy(this->orbitals[nHistory]);
    dPhi = orbital::param_copy(this->dOrbitals[nHistory]);

    int m = 0;
    for (int n = 0; n < nOrbitals; n++) {
        if (this->sepOrbitals) m = n;
        if (mpi::my_orb(Phi[n])) {
            std::vector<ComplexDouble> totCoefs;
            QMFunctionVector totOrbs;

            Orbital &phi_m = this->orbitals[nHistory][n];
            Orbital &fPhi_m = this->dOrbitals[nHistory][n];
            totCoefs.push_back(1.0);
            totOrbs.push_back(fPhi_m);

            // Ref. Harrisons KAIN paper the following has the wrong sign,
            // but we define the updates (lowercase f) with opposite sign
            // (but not the orbitals themselves).
            for (int j = 0; j < nHistory; j++) {
                ComplexVector partCoefs(4);
                QMFunctionVector partOrbs;

                partCoefs(0) = 1.0;
                Orbital &phi_j = this->orbitals[j][n];
                if (not mpi::my_orb(phi_j)) MSG_FATAL("MPI rank mismatch: phi_j");
                partOrbs.push_back(phi_j);

                partCoefs(1) = 1.0;
                Orbital &fPhi_j = this->dOrbitals[j][n];
                if (not mpi::my_orb(fPhi_j)) MSG_FATAL("MPI rank mismatch: fPhi_j");
                partOrbs.push_back(fPhi_j);

                partCoefs(2) = -1.0;
                partOrbs.push_back(phi_m);

                partCoefs(3) = -1.0;
                partOrbs.push_back(fPhi_m);

                Orbital partStep = phi_m.paramCopy();
                qmfunction::linear_combination(partStep, partCoefs, partOrbs, prec);

                ComplexDouble c_j = this->c[m](j);
                totCoefs.push_back(c_j);
                totOrbs.push_back(partStep);
            }

            // std::vector -> ComplexVector
            ComplexVector coefsVec(totCoefs.size());
            for (int i = 0; i < totCoefs.size(); i++) coefsVec(i) = totCoefs[i];

            dPhi[n] = Phi[n].paramCopy();
            qmfunction::linear_combination(dPhi[n], coefsVec, totOrbs, prec);
        }
    }

    // Treat Fock matrix as a whole using Frobenius inner product
    if (this->fock.size() == this->orbitals.size()) {
        if (F == 0 or dF == 0) MSG_ERROR("Invalid fock matrix");

        ComplexMatrix X_m = this->fock[nHistory];
        const ComplexMatrix &fX_m = this->dFock[nHistory];
        ComplexMatrix fockStep = ComplexMatrix::Zero(nOrbitals, nOrbitals);
        fockStep = fX_m;
        m = this->c.size();
        for (int j = 0; j < nHistory; j++) {
            const ComplexMatrix &X_j = this->fock[j];
            const ComplexMatrix &fX_j = this->dFock[j];
            ComplexMatrix tmpX = ComplexMatrix::Zero(nOrbitals, nOrbitals);
            tmpX += X_j;
            tmpX -= X_m;
            tmpX += fX_j;
            tmpX -= fX_m;
            tmpX *= this->c[m - 1](j);
            fockStep += tmpX;
        }
        *F = X_m;
        *dF = fockStep;
    }
    timer.stop();
    double t = timer.getWallTime();
    Printer::printDouble(0, "Expand solution", t, 5);
}

} // namespace mrchem
