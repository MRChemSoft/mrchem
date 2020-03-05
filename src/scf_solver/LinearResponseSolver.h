/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "SCFSolver.h"

/** @class LinearResponseSolver
 *
 */

namespace mrchem {

class Molecule;
class FockOperator;

class LinearResponseSolver final : public SCFSolver {
public:
    explicit LinearResponseSolver(bool dyn = false)
            : dynamic(dyn) {}
    ~LinearResponseSolver() override = default;

    bool optimize(double omega, Molecule &mol, FockOperator &F_0, FockOperator &F_1);
    void setOrthPrec(double prec) { this->orth_prec = prec; }

protected:
    const bool dynamic;
    double orth_prec{mrcpp::MachineZero};

    void printProperty() const;
    void printParameters(double omega, const std::string &oper) const;
};

} // namespace mrchem
