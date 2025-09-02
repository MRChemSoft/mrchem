
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

#include "qmoperators/QMPotential.h"
#include "tensor/RankOneOperator.h"
#include "tensor/RankZeroOperator.h"
#include <string>
#include "chemistry/Nucleus.h"
#include "pseudopotential/projectorOperator.h"

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
class AZoraPotential;
class NuclearOperator;
class CoulombOperator;
class ExchangeOperator;
class XCOperator;
class ElectricFieldOperator;
class ReactionOperator;

class FockBuilder final {
public:
    MomentumOperator &momentum() { return *this->mom; }
    RankZeroOperator &potential() { return this->V; }
    RankZeroOperator &perturbation() { return this->H_1; }

    std::shared_ptr<MomentumOperator> &getMomentumOperator() { return this->mom; }
    std::shared_ptr<NuclearOperator> &getNuclearOperator() { return this->nuc; }
    std::shared_ptr<CoulombOperator> &getCoulombOperator() { return this->coul; }
    std::shared_ptr<ExchangeOperator> &getExchangeOperator() { return this->ex; }
    std::shared_ptr<XCOperator> &getXCOperator() { return this->xc; }
    std::shared_ptr<ElectricFieldOperator> &getExtOperator() { return this->ext; }
    std::shared_ptr<ReactionOperator> &getReactionOperator() { return this->Ro; }
    std::shared_ptr<AZoraPotential> &getAZoraChiPotential() { return this->chiPot; }
    std::shared_ptr<ProjectorOperator> &getProjectorOperator() { return this->pp_projector; }

    void rotate(const ComplexMatrix &U);

    void build(double exx = 1.0);
    void setup(double prec);
    void clear();

    void setLightSpeed(double c) { this->light_speed = c; }
    double getLightSpeed() const { return this->light_speed; }

    bool isAZora() const { return zora_is_azora; }
    bool isZora() const { return (zora_has_nuc || zora_has_coul || zora_has_xc); }
    void setZoraType(bool has_nuc, bool has_coul, bool has_xc, bool is_azora);
    void setAZORADirectory(const std::string &dir) {azora_dir = dir;}
    void setNucs(const Nuclei &nucs) { this->nucs = nucs;}

    SCFEnergy trace(OrbitalVector &Phi, const Nuclei &nucs);
    ComplexMatrix operator()(OrbitalVector &bra, OrbitalVector &ket);
    ComplexMatrix kineticMatrix(OrbitalVector &bra, OrbitalVector &ket);
    ComplexMatrix potentialMatrix(OrbitalVector &bra, OrbitalVector &ket);

    OrbitalVector buildHelmholtzArgument(double prec, OrbitalVector Phi, ComplexMatrix F_mat, ComplexMatrix L_mat);

private:
    bool zora_has_nuc{false};
    bool zora_has_coul{false};
    bool zora_has_xc{false};
    bool zora_is_azora{false};
    std::string azora_dir = "";
    std::string azora_dir_src = "";
    std::string azora_dir_install = "";

    double light_speed{-1.0};
    double exact_exchange{1.0};
    RankZeroOperator zora_base;

    double prec;
    Nuclei nucs;

    RankZeroOperator V;   ///< Total potential energy operator
    RankZeroOperator H_1; ///< Perturbation operators

    std::shared_ptr<MomentumOperator> mom{nullptr};
    std::shared_ptr<NuclearOperator> nuc{nullptr};
    std::shared_ptr<CoulombOperator> coul{nullptr};
    std::shared_ptr<ExchangeOperator> ex{nullptr};
    std::shared_ptr<XCOperator> xc{nullptr};
    std::shared_ptr<ReactionOperator> Ro{nullptr};       // Reaction field operator
    std::shared_ptr<ElectricFieldOperator> ext{nullptr}; // Total external potential
    std::shared_ptr<ZoraOperator> chi{nullptr};
    std::shared_ptr<ZoraOperator> chi_inv{nullptr};
    std::shared_ptr<ProjectorOperator> pp_projector{nullptr};


    std::shared_ptr<QMPotential> collectZoraBasePotential();
    OrbitalVector buildHelmholtzArgumentZORA(OrbitalVector &Phi, OrbitalVector &Psi, DoubleVector eps, double prec);
    OrbitalVector buildHelmholtzArgumentNREL(OrbitalVector &Phi, OrbitalVector &Psi);
    std::shared_ptr<AZoraPotential> chiPot{nullptr}; // Potential for AZORA chi operator
    std::shared_ptr<QMPotential> chiInvPot{nullptr}; // Potential for AZORA chi_inv operator
};

} // namespace mrchem
