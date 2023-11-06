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

#include "ReactionPotential.h"
#include "environment/SCRF.h"

namespace mrchem {
class ReactionPotentialD2 final : public ReactionPotential {
public:
    ReactionPotentialD2(std::unique_ptr<SCRF> scrf, std::shared_ptr<mrchem::OrbitalVector> Phi, std::shared_ptr<OrbitalVector> X, std::shared_ptr<OrbitalVector> Y, bool mpi_share = false)
            : ReactionPotential(std::move(scrf), Phi, mpi_share)
            , orbitals_x(X)
            , orbitals_y(Y) {}

private:
    std::shared_ptr<OrbitalVector> orbitals_x; ///< Perturbed orbitals
    std::shared_ptr<OrbitalVector> orbitals_y; ///< Perturbed orbitals

    mrcpp::ComplexFunction &computePotential(double prec) const override;
};
} // namespace mrchem
