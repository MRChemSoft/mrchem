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

#pragma once

#include <string>
#include <vector>

#include "mrchem.h"
#include "qmfunctions/qmfunction_fwd.h"

#include "qmoperators/one_electron/ZoraOperator.h"

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

class SCFSolver {
public:
    SCFSolver() = default;
    virtual ~SCFSolver() = default;

    void setHistory(int hist) { this->history = hist; }
    void setCheckpoint(bool chk) { this->checkpoint = chk; }
    void setThreshold(double orb, double prop);
    void setOrbitalPrec(double init, double final);
    void setHelmholtzPrec(double prec) { this->helmPrec = prec; }
    void setMaxIterations(int iter) { this->maxIter = iter; }
    void setMethodName(const std::string &name) { this->methodName = name; }
    void setRelativityName(const std::string &name) { this->relativityName = name; }
    void setEnvironmentName(const std::string &name) { this->environmentName = name; }
    void setExternalFieldName(const std::string &name) { this->externalFieldName = name; }

protected:
    int history{0};                        ///< Maximum length of KAIN history
    int maxIter{-1};                       ///< Maximum number of iterations
    bool checkpoint{false};                ///< Dump orbitals to file every iteration
    double orbThrs{-1.0};                  ///< Convergence threshold for norm of orbital update
    double propThrs{-1.0};                 ///< Convergence threshold for property
    double helmPrec{-1.0};                 ///< Precision for construction of Helmholtz operators
    double orbPrec[3]{-1.0, -1.0, -1.0};   ///< Dynamic precision: [current_prec, start_prec, end_prec]
    std::string methodName;                ///< Name of electronic structure method to appear in output
    std::string relativityName{"None"};    ///< Name of ZORA method
    std::string environmentName{"None"};   ///< Name for external environment
    std::string externalFieldName{"None"}; ///< Name for external fields

    std::vector<double> error;    ///< Convergence orbital error
    std::vector<double> property; ///< Convergence property error

    virtual void reset();

    double adjustPrecision(double error);
    double getHelmholtzPrec();

    double getUpdate(const std::vector<double> &vec, int i, bool absPrec) const;
    void printUpdate(int plevel, const std::string &txt, double P, double dP, double thrs) const;

    bool checkConvergence(double err_o, double err_p) const;
    void printConvergence(bool converged, const std::string &txt) const;
    void printConvergenceHeader(const std::string &txt) const;
    void printConvergenceRow(int i) const;
    void printOrbitals(const DoubleVector &norms, const DoubleVector &errors, OrbitalVector &Phi, int flag, bool print_head = true) const;
    void printResidual(double residual, bool converged) const;
    void printMemory() const;
};

} // namespace mrchem
