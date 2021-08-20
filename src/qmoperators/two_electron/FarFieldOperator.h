/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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
#include "FarFieldPotential.h"
#include "chemistry/chemistry_utils.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/qmfunction_utils.h"

#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"

/** @class getFarFieldOperator
 *
 * @brief Operator containing a single CoulombPotential
 *
 * This class is a simple TensorOperator realization of @class CoulombPotential.
 *
 */

namespace mrchem {

class FarFieldOperator final : public RankZeroOperator {
public:
    FarFieldOperator(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi, const Nuclei &nucs, double exp_prec, bool mpi_share = false) {
        potential = std::make_shared<FarFieldPotential>(P, Phi, nucs, exp_prec, mpi_share);

        // Invoke operator= to assign *this operator
        RankZeroOperator &J = (*this);
        J = potential;
    }

    ~FarFieldOperator() override = default;

    auto &getPoisson() { return this->potential->getPoisson(); }
    auto &getDensity() { return this->potential->getDensity(); }
    const Nuclei &getNuclei() const { return this->potential->getNuclei(); }

    double getNucPrec() { return this->potential->getNucPrec(); }

    ComplexDouble trace(OrbitalVector &Phi) { return 0.5 * RankZeroOperator::trace(Phi); }
    using RankZeroOperator::trace;

private:
    std::shared_ptr<FarFieldPotential> potential{nullptr};
};

} // namespace mrchem
