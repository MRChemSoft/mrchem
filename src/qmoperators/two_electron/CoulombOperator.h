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

#include "tensor/RankZeroOperator.h"

#include "CoulombPotential.h"
#include "CoulombPotentialD1.h"
#include "CoulombPotentialD2.h"

/** @class CoulombOperator
 *
 * @brief Operator containing a single CoulombPotential
 *
 * This class is a simple TensorOperator realization of @class CoulombPotential.
 *
 */

namespace mrchem {

class CoulombOperator final : public RankZeroOperator {
public:
    explicit CoulombOperator(std::shared_ptr<mrcpp::PoissonOperator> P, bool mpi_share = false) {
        potential = std::make_shared<CoulombPotential>(P, nullptr, mpi_share);

        // Invoke operator= to assign *this operator
        RankZeroOperator &J = (*this);
        J = potential;
        J.name() = "J";
    }
    CoulombOperator(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi, bool mpi_share = false) {
        potential = std::make_shared<CoulombPotentialD1>(P, Phi, mpi_share);

        // Invoke operator= to assign *this operator
        RankZeroOperator &J = (*this);
        J = potential;
        J.name() = "J";
    }
    CoulombOperator(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi, std::shared_ptr<OrbitalVector> X, std::shared_ptr<OrbitalVector> Y, bool mpi_share = false) {
        potential = std::make_shared<CoulombPotentialD2>(P, Phi, X, Y, mpi_share);

        // Invoke operator= to assign *this operator
        RankZeroOperator &J = (*this);
        J = potential;
        J.name() = "J";
    }
    ~CoulombOperator() override = default;

    auto &getPoisson() { return this->potential->getPoisson(); }
    auto &getDensity() { return this->potential->getDensity(); }

private:
    std::shared_ptr<CoulombPotential> potential{nullptr};
};

} // namespace mrchem
