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

#include "ReactionPotential.h"
#include "qmoperators/RankZeroTensorOperator.h"

/** @class ReactionOperator
 *
 * @brief Operator containing a single ReactionPotential
 *
 * This class is a simple TensorOperator realization of @class ReactionPotential.
 *
 */

namespace mrchem {
class SCRF;

class ReactionOperator final : public RankZeroTensorOperator {
public:
    ReactionOperator(std::shared_ptr<mrchem::OrbitalVector> Phi_p, SCRF help) {
        potential = std::make_shared<ReactionPotential>(Phi_p, help);
        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
    }

    ComplexDouble trace(OrbitalVector &Phi) { return RankZeroTensorOperator::trace(Phi); }

    double getTotalEnergy() { return this->potential->getTotalEnergy(); }
    double getElectronicEnergy() { return this->potential->getElectronicEnergy(); }
    double getNuclearEnergy() { return this->potential->getNuclearEnergy(); }
    SCRF getHelper() { return this->potential->getHelper(); }
    std::shared_ptr<ReactionPotential> getPotential() { return this->potential; }
    void updateMOResidual(double const err_t) { this->potential->updateMOResidual(err_t); }

    QMFunction &getCurrentReactionPotential() { return this->potential->getCurrentReactionPotential(); }
    QMFunction &getPreviousReactionPotential() { return this->potential->getPreviousReactionPotential(); }
    QMFunction &getCurrentDifferenceReactionPotential() {
        return this->potential->getCurrentDifferenceReactionPotential();
    }

    QMFunction &getCurrentGamma() { return this->potential->getCurrentGamma(); }
    QMFunction &getPreviousGamma() { return this->potential->getPreviousGamma(); }
    QMFunction &getCurrentDifferenceGamma() { return this->potential->getCurrentDifferenceGamma(); }

    void setTesting() { this->potential->setTesting(); }

private:
    std::shared_ptr<ReactionPotential> potential{nullptr};
};

} // namespace mrchem
