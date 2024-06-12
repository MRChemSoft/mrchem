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

#include "XCPotential.h"
#include "XCPotentialD1.h"
#include "XCPotentialD2.h"

/** @class XCOperator
 *
 * @brief DFT Exchange-Correlation operator containing a single XCPotential
 *
 * This class is a simple TensorOperator realization of @class XCPotential.
 *
 */

namespace mrchem {

class XCOperator final : public RankZeroOperator {
public:
    explicit XCOperator(std::unique_ptr<mrdft::MRDFT> &F, std::shared_ptr<OrbitalVector> Phi = nullptr, bool mpi_shared = false) {
        potential = std::make_shared<XCPotentialD1>(F, Phi, mpi_shared);

        // Invoke operator= to assign *this operator
        RankZeroOperator &XC = (*this);
        XC = potential;
        XC.name() = "V_xc";
    }
    XCOperator(std::unique_ptr<mrdft::MRDFT> &F, std::shared_ptr<OrbitalVector> Phi, std::shared_ptr<OrbitalVector> X, std::shared_ptr<OrbitalVector> Y, bool mpi_shared = false) {
        potential = std::make_shared<XCPotentialD2>(F, Phi, X, Y, mpi_shared);

        // Invoke operator= to assign *this operator
        RankZeroOperator &XC = (*this);
        XC = potential;
        XC.name() = "V_xc";
    }
    ~XCOperator() override = default;

    auto getEnergy() { return potential->getEnergy(); }
    auto &getDensity(DensityType spin, int pert_idx = 0) { return potential->getDensity(spin, pert_idx); }

    void setSpin(int spin) {
        mrcpp::FunctionTree<3> &pot = this->potential->getPotential(spin);
        this->potential->setReal(&pot);
    }
    void clearSpin() { this->potential->setReal(nullptr); }

    std::shared_ptr<XCPotential> potential{nullptr};
};

} // namespace mrchem
