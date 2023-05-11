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

#include <memory>
#include <MRCPP/Parallel>

#include "qmoperators/QMOperator.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/qmfunction_fwd.h"

namespace mrchem {

/** @class ExchangePotential
 *
 *  @brief Hartree-Fock exchange potential defined by a particular set of orbitals
 *
 * The operator is defined as the Hartree-Fock exchange arising from a particular
 * set of orbitals. The OrbitalVector defining the operator is fixed throughout the
 * operators life time, but the orbitals themselves are allowed to change in between
 * each application. The internal exchange potentials (the operator applied to it's
 * own orbitals) can be precomputed and stored for fast retrieval. Option to use
 * screening based on previous calculations of the internal exchange (make sure that
 * the internal orbitals haven't been significantly changed since the last time the
 * operator was set up, e.g. through an orbital rotation).
 */

class ExchangePotential : public QMOperator {
public:
    ExchangePotential(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi, double prec);
    ~ExchangePotential() override = default;

    friend class ExchangeOperator;

protected:
    bool pre_compute{false};                         ///< Precompute internal exchange
    double exchange_prec;                            ///< Screening precision for exchange construction
    OrbitalVector exchange;                          ///< Precomputed exchange from the internal orbital set
    std::shared_ptr<OrbitalVector> orbitals;         ///< Internal orbitals defining the exchange operator
    std::shared_ptr<mrcpp::PoissonOperator> poisson; ///< Poisson operator to compute orbital contributions

    void setPreCompute() { this->pre_compute = true; }

    auto &getPoisson() { return this->poisson; }
    double getSpinFactor(Orbital phi_i, Orbital phi_j) const;

    void rotate(const ComplexMatrix &U);
    void setup(double prec) override;
    void clear() override;

    virtual void setupBank() = 0;
    virtual void clearBank() {}

    virtual int testInternal(Orbital phi_p) const { return -1; }
    virtual void setupInternal(double prec) {}
    void clearInternal() { this->exchange.clear(); }

    void calcExchange_kij(double prec, Orbital phi_k, Orbital phi_i, Orbital phi_j, Orbital &out_kij, Orbital *out_jji = nullptr);
};

} // namespace mrchem
