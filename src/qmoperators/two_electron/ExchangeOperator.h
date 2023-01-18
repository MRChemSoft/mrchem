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

#include "ExchangePotentialD1.h"
#include "ExchangePotentialD2.h"

/**
 * @class ExchangeOperator
 *
 * @brief Operator containing a single ExchangePotential
 *
 * The operator is defined as the Hartree-Fock exchange arising from a particular
 * set of orbitals. The OrbitalVector defining the operator is fixed throughout the
 * operators life time, but the orbitals themselves are allowed to change in between
 * each application. The internal exchange potentials (the operator applied to it's
 * own orbitals) can be precomputed and stored for fast retrieval.
 */

namespace mrchem {

class ExchangeOperator final : public RankZeroOperator {
public:
    ExchangeOperator(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi, double exchange_prec = -1.0) {
        exchange = std::make_shared<ExchangePotentialD1>(P, Phi, exchange_prec);

        // Invoke operator= to assign *this operator
        RankZeroOperator &K = (*this);
        K = exchange;
        K.name() = "K";
    }

    ExchangeOperator(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi, std::shared_ptr<OrbitalVector> X, std::shared_ptr<OrbitalVector> Y, double exchange_prec = -1.0) {
        exchange = std::make_shared<ExchangePotentialD2>(P, Phi, X, Y, exchange_prec);

        // Invoke operator= to assign *this operator
        RankZeroOperator &K = (*this);
        K = exchange;
        K.name() = "K";
    }

    ~ExchangeOperator() override = default;

    auto &getPoisson() { return exchange->getPoisson(); }
    void setPreCompute() { exchange->setPreCompute(); }
    void rotate(const ComplexMatrix &U) { exchange->rotate(U); }

    ComplexDouble trace(OrbitalVector &Phi) { return 0.5 * RankZeroOperator::trace(Phi); }

private:
    std::shared_ptr<ExchangePotential> exchange{nullptr};
};

} // namespace mrchem
