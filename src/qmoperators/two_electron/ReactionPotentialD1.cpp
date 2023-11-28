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

#include "ReactionPotentialD1.h"

#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "qmfunctions/density_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator = mrcpp::PoissonOperator;
using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;

namespace mrchem {
mrcpp::ComplexFunction &ReactionPotentialD1::computePotential(double prec) const {
    // construct electronic density from the orbitals
    OrbitalVector &Phi = *this->orbitals;
    Density rho_el(false);
    density::compute(this->apply_prec, rho_el, Phi, DensityType::Total);
    // change sign, because it's the electronic density
    rho_el.rescale(-1.0);

    return this->solver->solveEquation(prec, rho_el);
}
} // namespace mrchem
