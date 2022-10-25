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

#include "HartreePotential.h"
#include "tensor/RankZeroOperator.h"

/** @class CoulombOperator
 *
 * @brief Operator containing a single CoulombPotential
 *
 * This class is a simple TensorOperator realization of @class CoulombPotential.
 *
 */

namespace mrchem {

class HartreeOperator final : public RankZeroOperator {
public:
    HartreeOperator(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi, const Nuclei &nucs, const double &rc) {

        potential = std::make_shared<HartreePotential>(P, Phi, nucs, rc);
        RankZeroOperator &J = (*this);
        J = potential;
    }

    ~HartreeOperator() override = default;

    auto &getPoisson() { return this->potential->getPoisson(); }
    auto &getDensity() { return this->potential->getDensity(); }
    auto &getPotential() { return this->potential; }

private:
    std::shared_ptr<HartreePotential> potential{nullptr};
};

} // namespace mrchem
