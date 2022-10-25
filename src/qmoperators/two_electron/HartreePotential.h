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

#include "chemistry/Nucleus.h"
#include "qmfunctions/Density.h"
#include "qmoperators/QMPotential.h"

namespace mrchem {

class HartreePotential final : public QMPotential {
public:
    HartreePotential(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi, const Nuclei &nucs, const double &rc);
    ~HartreePotential() override = default;

    friend class HartreeOperator;

    Density &getBSmear() { return this->b_smeared; }
    double getRc() { return this->rc; }

private:
    Density density;                                 ///< Ground-state electron density
    std::shared_ptr<OrbitalVector> orbitals;         ///< Unperturbed orbitals defining the ground-state electron density
    std::shared_ptr<mrcpp::PoissonOperator> poisson; ///< Operator used to compute the potential
    Density b_smeared{false};
    Nuclei nuclei;
    double rc;

    auto &getPoisson() { return this->poisson; }
    auto &getDensity() { return this->density; }

    bool hasDensity() const { return (this->density.squaredNorm() < 0.0) ? false : true; }

    void setup(double prec) override;
    void clear() override;

    void setupDensity(double prec);
    void setupPotential(double prec);
};

} // namespace mrchem
