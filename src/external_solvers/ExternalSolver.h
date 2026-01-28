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

#include <vector>

#include "mrchem.h"
#include "qmfunctions/Orbital.h"
#include "qmoperators/two_electron/FockBuilder.h"

// using PoissonOperator = mrcpp::PoissonOperator;

namespace mrchem {

class FockBuilder;

class ExternalSolver {
public:
    ExternalSolver() = default;
    virtual ~ExternalSolver() = default;

    double get_energy() { return this->energy; }
    // void set_precision(double p){this->prec=p;}

    void set_integrals(OrbitalVector &Phi, FockBuilder &F);
    std::shared_ptr<ComplexMatrix> get_one_body_integrals() { return this->one_body_integrals; }
    std::shared_ptr<ComplexTensorR4> get_two_body_integrals() { return this->two_body_integrals; }
    std::shared_ptr<ComplexMatrix> get_one_rdm() { return this->one_rdm; }
    std::shared_ptr<ComplexTensorR4> get_two_rdm() { return this->two_rdm; }

    virtual void optimize() = 0;

protected:
    // Dummy default precision
    double prec = 1e-3;
    double energy{};

    std::shared_ptr<ComplexMatrix> one_body_integrals{};
    std::shared_ptr<ComplexTensorR4> two_body_integrals{};
    std::shared_ptr<ComplexMatrix> one_rdm{};
    std::shared_ptr<ComplexTensorR4> two_rdm{};

    void set_one_body_integrals(OrbitalVector &Phi, MomentumOperator &P, NuclearOperator &V);
    void set_two_body_integrals(OrbitalVector &Phi, GenericTwoOrbitalsOperator &g);

    virtual void calculate_rdms() = 0;
};

} // namespace mrchem