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

/** @class FockOperator
 *
 * @brief Operator containing the standard SCF operators
 *
 * This is a simple collection of operators used in ground state SCF calculations.
 * The operator is separated into kinetic and potential parts, since the MW way of
 * solving the SCF equations is to invert the kinetic part, and apply the potential
 * part as usual.
 */

namespace mrchem {

class SCFEnergy;
class KineticOperator;
class NuclearOperator;
class SmearedNuclearOperator;
class CoulombOperator;
class FarFieldOperator;
class ExchangeOperator;
class XCOperator;
class ElectricFieldOperator;
class ReactionOperator;

class FockOperator final : public RankZeroOperator {
public:
    FockOperator(std::shared_ptr<KineticOperator> t = nullptr,
                 std::shared_ptr<NuclearOperator> v = nullptr,
                 std::shared_ptr<SmearedNuclearOperator> sv = nullptr,
                 std::shared_ptr<CoulombOperator> j = nullptr,
                 std::shared_ptr<FarFieldOperator> f = nullptr,
                 std::shared_ptr<ExchangeOperator> k = nullptr,
                 std::shared_ptr<XCOperator> xc = nullptr,
                 std::shared_ptr<ElectricFieldOperator> ext = nullptr,
                 std::shared_ptr<ReactionOperator> reo = nullptr);

    RankZeroOperator &kinetic() { return this->T; }
    RankZeroOperator &potential() { return this->V; }
    RankZeroOperator &perturbation() { return this->H_1; }

    std::shared_ptr<KineticOperator> &getKineticOperator() { return this->kin; }
    std::shared_ptr<NuclearOperator> &getNuclearOperator() { return this->nuc; }
    std::shared_ptr<SmearedNuclearOperator> &getSmearedNuclearOperator() { return this->snuc; }
    std::shared_ptr<CoulombOperator> &getCoulombOperator() { return this->coul; }
    std::shared_ptr<FarFieldOperator> &getFarFieldOperator() { return this->ff; }
    std::shared_ptr<ExchangeOperator> &getExchangeOperator() { return this->ex; }
    std::shared_ptr<XCOperator> &getXCOperator() { return this->xc; }
    std::shared_ptr<ElectricFieldOperator> &getExtOperator() { return this->ext; }
    std::shared_ptr<ReactionOperator> &getReactionOperator() { return this->Ro; }

    void rotate(const ComplexMatrix &U);

    void build(double exx = 1.0);
    void setup(double prec);
    void clear();

    SCFEnergy trace(OrbitalVector &Phi, const Nuclei &nucs);

    ComplexMatrix operator()(OrbitalVector &bra, OrbitalVector &ket);
    ComplexMatrix dagger(OrbitalVector &bra, OrbitalVector &ket);

    using RankZeroOperator::operator();
    using RankZeroOperator::dagger;

private:
    double exact_exchange{1.0};
    RankZeroOperator T;   ///< Total kinetic energy operator
    RankZeroOperator V;   ///< Total potential energy operator
    RankZeroOperator H_1; ///< Perturbation operators

    std::shared_ptr<KineticOperator> kin;
    std::shared_ptr<NuclearOperator> nuc;
    std::shared_ptr<SmearedNuclearOperator> snuc;
    std::shared_ptr<CoulombOperator> coul;
    std::shared_ptr<FarFieldOperator> ff;
    std::shared_ptr<ExchangeOperator> ex;
    std::shared_ptr<XCOperator> xc;
    std::shared_ptr<ElectricFieldOperator> ext; ///< Total external potential
    std::shared_ptr<ReactionOperator> Ro;       ///< Reaction field operator
};

} // namespace mrchem
