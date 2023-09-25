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

#include "ReactionPotentialD1.h"
#include "ReactionPotentialD2.h"
#include "environment/GPESolver.h"

/** @class ReactionOperator
 *
 * @brief Operator containing a single ReactionPotential
 *
 * This class is a simple TensorOperator realization of @class ReactionPotential.
 *
 */

namespace mrchem {
class ReactionOperator final : public RankZeroOperator {
public:
    ReactionOperator(std::unique_ptr<GPESolver> gpesolver_p, std::shared_ptr<mrchem::OrbitalVector> Phi_p, bool mpi_share = false) {
        potential = std::make_shared<ReactionPotentialD1>(std::move(gpesolver_p), Phi_p, mpi_share);
        // Invoke operator= to assign *this operator
        RankZeroOperator &V = (*this);
        V = potential;
        V.name() = "V_r";
    }

    ReactionOperator(std::unique_ptr<GPESolver> gpesolver_p, std::shared_ptr<mrchem::OrbitalVector> Phi, std::shared_ptr<OrbitalVector> X, std::shared_ptr<OrbitalVector> Y, bool mpi_share = false) {
        // check that the GPESolver object uses the electronic density only
        if (gpesolver_p->getDensityType() != SCRFDensityType::ELECTRONIC) MSG_ERROR("Invalid SCRF object passed: only electronic density in response");
        potential = std::make_shared<ReactionPotentialD2>(std::move(gpesolver_p), Phi, X, Y, mpi_share);
        // Invoke operator= to assign *this operator
        RankZeroOperator &V = (*this);
        V = potential;
        V.name() = "V_r";
    }

    ComplexDouble trace(OrbitalVector &Phi) { return RankZeroOperator::trace(Phi); }

    GPESolver *getSolver() { return this->potential->getSolver(); }
    std::shared_ptr<ReactionPotential> getPotential() { return this->potential; }
    void updateMOResidual(double const err_t) { this->potential->updateMOResidual(err_t); }

private:
    std::shared_ptr<ReactionPotential> potential{nullptr};
};
} // namespace mrchem
