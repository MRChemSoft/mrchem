/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include <nlohmann/json.hpp>

#include "SCFSolver.h"

/** @class ExcitedStatesSolver
 *
 */

namespace mrchem {

class Molecule;
class FockBuilder;

class ExcitedStatesSolver final : public SCFSolver {
public:
    explicit ExcitedStatesSolver(bool dyn = false, int N_states = 1)
            : dynamic(dyn)
            , n_states(N_states) {}
    ~ExcitedStatesSolver() override = default;

    nlohmann::json optimize(Molecule &mol, FockBuilder &F_0, std::vector<FockBuilder> &F_1_vec);
    void setOrthPrec(double prec) { this->orth_prec = prec; }
    void setCheckpointFile(const std::string &file_x, const std::string &file_y) {
        this->chkFileX = file_x;
        this->chkFileY = file_y;
    }

protected:
    const bool dynamic;
    const int n_states;
    double orth_prec{mrcpp::MachineZero};
    std::string chkFileX; ///< Name of checkpoint file
    std::string chkFileY; ///< Name of checkpoint file

    double computeOmega(OrbitalVector &Phi_0, OrbitalVector &X, FockBuilder &F_0, RankZeroOperator &V_1, ComplexMatrix &F_mat_0);
    double updateOmega(OrbitalVector &X_n, OrbitalVector &X_np1);
    void printProperty() const;
    void printParameters(double omega, const std::string &oper) const;
};

} // namespace mrchem
