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

#pragma once

#include <string>
#include <vector>

#include "mrchem.h"
#include "qmfunctions/qmfunction_fwd.h"

/** @class SCF
 *
 * @brief Abstract base class for different types of SCF solvers
 *
 * The ground state and linear response SCF solvers share some common features which
 * are collected in this abstract base class. This mainly involves convergence routines
 * and printing.
 *
 */

namespace mrchem {

class HelmholtzVector;
class Accelerator;
class FockOperator;

class SCF {
public:
    SCF(HelmholtzVector &h);

    virtual bool optimize() = 0;

    void setRotation(int iter) { this->rotation = iter; }
    void setCanonical(bool can) { this->canonical = can; }
    void setThreshold(double orb, double prop);
    void setOrbitalPrec(double init, double final);
    void setMaxIterations(int m_iter) { this->maxIter = m_iter; }

protected:
    int maxIter;       ///< Maximum number of iterations
    int rotation;      ///< Number of iterations between localization/diagonalization
    bool canonical;    ///< Use localized or canonical orbitals
    double orbThrs;    ///< Convergence threshold for norm of orbital update
    double propThrs;   ///< Convergence threshold for property
    double orbPrec[3]; ///< Dynamic precision: [current_prec, start_prec, end_prec]

    std::vector<double> orbError; ///< Convergence orbital error
    std::vector<double> property; ///< Convergence property error

    HelmholtzVector *helmholtz; ///< Pointer to external object

    bool checkConvergence(double err_o, double err_p) const;
    bool needLocalization(int nIter) const;
    bool needDiagonalization(int nIter) const;

    double adjustPrecision(double error);
    void resetPrecision() { this->orbPrec[0] = this->orbPrec[1]; }

    double getUpdate(const std::vector<double> &vec, int i, bool absPrec) const;
    void printUpdate(const std::string &name, double P, double dP) const;

    void printOrbitals(const DoubleVector &epsilon, const OrbitalVector &Phi, int flag) const;
    void printConvergence(bool converged) const;
    void printCycleHeader(int nIter) const;
    void printCycleFooter(double t) const;
};

} // namespace mrchem
