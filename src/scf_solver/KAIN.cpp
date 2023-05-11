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

#include <MRCPP/Parallel>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "KAIN.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

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
    Timer t_tot;
    int nHistory = this->orbitals.size() - 1;
    if (nHistory < 1) MSG_ABORT("Not enough history to setup system of equations");

    std::vector<ComplexMatrix> A_matrices;
    std::vector<ComplexVector> b_vectors;

    int nOrbitals = this->orbitals[nHistory].size();
    for (int n = 0; n < nOrbitals; n++) {
        auto orbA = ComplexMatrix::Zero(nHistory, nHistory).eval();
        auto orbB = ComplexVector::Zero(nHistory).eval();

        auto &phi_m = this->orbitals[nHistory][n];
        auto &fPhi_m = this->dOrbitals[nHistory][n];
        if (mrcpp::mpi::my_orb(phi_m)) {

            for (int i = 0; i < nHistory; i++) {
                auto &phi_i = this->orbitals[i][n];
                auto dPhi_im = phi_m.paramCopy();
                mrcpp::cplxfunc::add(dPhi_im, 1.0, phi_i, -1.0, phi_m, -1.0);

                for (int j = 0; j < nHistory; j++) {
                    auto &fPhi_j = this->dOrbitals[j][n];
                    auto dfPhi_jm = fPhi_m.paramCopy();
                    mrcpp::cplxfunc::add(dfPhi_jm, 1.0, fPhi_j, -1.0, fPhi_m, -1.0);

                    // Ref. Harrisons KAIN paper the following has the wrong sign,
                    // but we define the updates (lowercase f) with opposite sign.
                    orbA(i, j) -= orbital::dot(dPhi_im, dfPhi_jm);
                }
                orbB(i) += orbital::dot(dPhi_im, fPhi_m);
            }
        }
        double alpha = (this->scaling.size() == nOrbitals) ? scaling[n] : 1.0;
        A_matrices.push_back(orbA * alpha);
        b_vectors.push_back(orbB);
    }

    for (int i = 0; i < nOrbitals; i++) {
        mrcpp::mpi::allreduce_matrix(A_matrices[i], mrcpp::mpi::comm_wrk);
        mrcpp::mpi::allreduce_vector(b_vectors[i], mrcpp::mpi::comm_wrk);
    }

    // Fock matrix is treated as a whole using the Frobenius inner product
    if (this->orbitals.size() == this->fock.size()) {
        const auto &X_m = this->fock[nHistory];
        const auto &fX_m = this->dFock[nHistory];

        auto fockA = ComplexMatrix::Zero(nHistory, nHistory).eval();
        auto fockB = ComplexVector::Zero(nHistory).eval();

        for (int i = 0; i < nHistory; i++) {
            const auto &X_i = this->fock[i];
            auto dX_im = X_i - X_m;
            for (int j = 0; j < nHistory; j++) {
                const auto &fX_j = this->dFock[j];
                auto dfX_jm = fX_j - fX_m;
                auto prod = dX_im.adjoint() * dfX_jm;
                fockA(i, j) -= prod.trace();
            }
            auto prod = dX_im.adjoint() * fX_m;
            fockB(i) += prod.trace();
        }
        A_matrices.push_back(fockA);
        b_vectors.push_back(fockB);
    }

    sortLinearSystem(A_matrices, b_vectors);
    mrcpp::print::time(this->pl + 2, "Setup linear system", t_tot);
}

/** @brief Compute the next step for orbitals and orbital updates
 *
 * The next step \f$ \delta x^n \f$ is constructed from the solution
 * solution \f$ c \f$ of the linear problem \f$ Ac = b \f$ as:
 *
 * \f$ \delta x^n = f(x^n) + \sum_{j=1}^m c_j[(x^j-x^n)+(f(x^j)-f(x^n))]\f$
 */
void KAIN::expandSolution(double prec, OrbitalVector &Phi, OrbitalVector &dPhi, ComplexMatrix *F, ComplexMatrix *dF) {
    Timer t_tot;
    int nHistory = this->orbitals.size() - 1;
    int nOrbitals = this->orbitals[nHistory].size();

    // Orbitals are unchanged, updates are overwritten
    dPhi = orbital::param_copy(this->dOrbitals[nHistory]);

    int m = 0;
    for (int n = 0; n < nOrbitals; n++) {
        if (this->sepOrbitals) m = n;
        if (mrcpp::mpi::my_orb(Phi[n])) {
            std::vector<ComplexDouble> totCoefs;
            std::vector<mrcpp::ComplexFunction> totOrbs;

            auto &phi_m = this->orbitals[nHistory][n];
            auto &fPhi_m = this->dOrbitals[nHistory][n];
            totCoefs.push_back({1.0, 0.0});
            totOrbs.push_back(fPhi_m);

            // Ref. Harrisons KAIN paper the following has the wrong sign,
            // but we define the updates (lowercase f) with opposite sign
            // (but not the orbitals themselves).
            for (int j = 0; j < nHistory; j++) {
                ComplexVector partCoefs(4);
                std::vector<mrcpp::ComplexFunction> partOrbs;

                partCoefs(0) = {1.0, 0.0};
                auto &phi_j = this->orbitals[j][n];
                partOrbs.push_back(phi_j);

                partCoefs(1) = {1.0, 0.0};
                auto &fPhi_j = this->dOrbitals[j][n];
                partOrbs.push_back(fPhi_j);

                partCoefs(2) = {-1.0, 0.0};
                partOrbs.push_back(phi_m);

                partCoefs(3) = {-1.0, 0.0};
                partOrbs.push_back(fPhi_m);

                auto partStep = phi_m.paramCopy();
                mrcpp::cplxfunc::linear_combination(partStep, partCoefs, partOrbs, prec);

                auto c_j = this->c[m](j);
                totCoefs.push_back(c_j);
                totOrbs.push_back(partStep);
            }

            // std::vector -> ComplexVector
            ComplexVector coefsVec(totCoefs.size());
            for (int i = 0; i < totCoefs.size(); i++) coefsVec(i) = totCoefs[i];

            dPhi[n] = Phi[n].paramCopy();
            mrcpp::cplxfunc::linear_combination(dPhi[n], coefsVec, totOrbs, prec);
        }
    }

    // Treat Fock matrix as a whole using Frobenius inner product
    if (this->fock.size() == this->orbitals.size()) {
        if (F == nullptr or dF == nullptr) MSG_ERROR("Invalid fock matrix");

        auto X_m = this->fock[nHistory];
        const auto &fX_m = this->dFock[nHistory];
        auto fockStep = ComplexMatrix::Zero(nOrbitals, nOrbitals).eval();
        fockStep = fX_m;
        m = this->c.size();
        for (int j = 0; j < nHistory; j++) {
            const auto &X_j = this->fock[j];
            const auto &fX_j = this->dFock[j];
            auto tmpX = ComplexMatrix::Zero(nOrbitals, nOrbitals).eval();
            tmpX += X_j;
            tmpX -= X_m;
            tmpX += fX_j;
            tmpX -= fX_m;
            tmpX *= this->c[m - 1](j);
            fockStep += tmpX;
        }
        // F is unchanged, dF is overwritten
        *dF = fockStep;
    }
    mrcpp::print::time(this->pl + 2, "Expand solution", t_tot);
}

} // namespace mrchem
