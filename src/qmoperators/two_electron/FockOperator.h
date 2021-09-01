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
#include "qmoperators/QMPotential.h"

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
class MomentumOperator;
class KineticOperator;
class ZoraKineticOperator;
class ZoraOperator;
class NuclearOperator;
class CoulombOperator;
class ExchangeOperator;
class XCOperator;
class ElectricFieldOperator;
class ReactionOperator;

class FockOperator final : public RankZeroOperator {
public:
    bool isZora() const { return (this->vz != nullptr); }
    ZoraOperator &zora() { return *this->vz; }
    MomentumOperator &momentum() { return *this->mom; }
    RankZeroOperator &kinetic() { return this->T; }
    RankZeroOperator &potential() { return this->V; }
    RankZeroOperator &perturbation() { return this->H_1; }
    RankZeroOperator &sqrt_zora_pot() { return this->sqrt_zora; }
    RankZeroOperator &mod_zora_pot() { return this->mod_zora; }

    std::shared_ptr<MomentumOperator> &getMomentumOperator() { return this->mom; }
    std::shared_ptr<NuclearOperator> &getNuclearOperator() { return this->nuc; }
    std::shared_ptr<CoulombOperator> &getCoulombOperator() { return this->coul; }
    std::shared_ptr<ExchangeOperator> &getExchangeOperator() { return this->ex; }
    std::shared_ptr<XCOperator> &getXCOperator() { return this->xc; }
    std::shared_ptr<ZoraOperator> &getZoraOperator() { return this->vz; }
    std::shared_ptr<ElectricFieldOperator> &getExtOperator() { return this->ext; }
    std::shared_ptr<ReactionOperator> &getReactionOperator() { return this->Ro; }
    std::shared_ptr<RankZeroOperator> &getSqrtZora() { return this->sqrt_vz; }
    std::shared_ptr<RankZeroOperator> &getModZora() { return this->mod_vz; }

    void rotate(const ComplexMatrix &U);

    void build(double exx = 1.0);
    void setup(double prec);
    void clear();

    SCFEnergy trace(OrbitalVector &Phi, const Nuclei &nucs);

    ComplexMatrix operator()(OrbitalVector &bra, OrbitalVector &ket);
    ComplexMatrix dagger(OrbitalVector &bra, OrbitalVector &ket);

    RankZeroOperator buildHelmholtzArgumentOperator(double prec);

    using RankZeroOperator::operator();
    using RankZeroOperator::dagger;

private:
    double exact_exchange{1.0};
    RankZeroOperator T;       ///< Total kinetic energy operator
    RankZeroOperator V;       ///< Total potential energy operator
    RankZeroOperator H_1;     ///< Perturbation operators
    RankZeroOperator sqrt_zora;
    RankZeroOperator mod_zora;

    std::shared_ptr<MomentumOperator> mom{nullptr};
    std::shared_ptr<ZoraOperator> vz{nullptr};
    std::shared_ptr<NuclearOperator> nuc{nullptr};
    std::shared_ptr<CoulombOperator> coul{nullptr};
    std::shared_ptr<ExchangeOperator> ex{nullptr};
    std::shared_ptr<XCOperator> xc{nullptr};
    std::shared_ptr<ReactionOperator> Ro{nullptr};       ///< Reaction field operator
    std::shared_ptr<ElectricFieldOperator> ext{nullptr}; ///< Total external potential
    std::shared_ptr<RankZeroOperator> sqrt_vz{nullptr};
    std::shared_ptr<RankZeroOperator> mod_vz{nullptr};
};

} // namespace mrchem
