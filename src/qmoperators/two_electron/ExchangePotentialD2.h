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

#include "ExchangePotential.h"

#include "qmfunctions/qmfunction_fwd.h"
#include "utils/Bank.h"

namespace mrchem {

/** @class ExchangePotentialD2
 *
 *  @brief Hartree-Fock exchange potential defined by a set of unperturbed orbitals
 *
 * The operator is defined as the Hartree-Fock exchange arising from a
 * set of unperturbed orbitals. The OrbitalVector defining the
 * operator is fixed throughout the operator life time, but the
 * orbitals themselves are allowed to change in between each
 * application. The internal exchange potentials (the operator applied
 * to it's own orbitals) can be precomputed and stored for fast
 * retrieval. Option to use screening based on previous calculations
 * of the internal exchange (make sure that the internal orbitals
 * haven't been significantly changed since the last time the operator
 * was set up, e.g. through an orbital rotation).
 */

class ExchangePotentialD2 final : public ExchangePotential {
public:
    ExchangePotentialD2(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi, std::shared_ptr<OrbitalVector> X, std::shared_ptr<OrbitalVector> Y, double prec);
    ~ExchangePotentialD2() override = default;

    friend class ExchangeOperator;

private:
    BankAccount PhiBank; // to put the Orbitals
    BankAccount XBank;
    BankAccount YBank;
    bool useOnlyX;                             ///< true if X and Y are the same set of orbitals
    std::shared_ptr<OrbitalVector> orbitals_x; ///< first set of perturbed orbitals defining the exchange operator
    std::shared_ptr<OrbitalVector> orbitals_y; ///< second set of perturbed orbitals defining the exchange operator

    void setupBank() override;
    void clearBank();

    ComplexDouble evalf(const mrcpp::Coord<3> &r) const override { return 0.0; }

    Orbital apply(Orbital phi_p) override;
    Orbital dagger(Orbital phi_p) override;
    QMOperatorVector apply(std::shared_ptr<QMOperator> &O) override;
};

} // namespace mrchem
