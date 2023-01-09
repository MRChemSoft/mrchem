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

#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombPotentialD2.h"
#include "qmfunctions/density_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator = mrcpp::PoissonOperator;
using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {

CoulombPotentialD2::CoulombPotentialD2(PoissonOperator_p P, OrbitalVector_p Phi, OrbitalVector_p X, OrbitalVector_p Y, bool mpi_share)
        : CoulombPotential(P, Phi, mpi_share)
        , orbitals_x(X)
        , orbitals_y(Y) {}

void CoulombPotentialD2::setupGlobalDensity(double prec) {
    if (hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_x == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_y == nullptr) MSG_ERROR("Orbitals not initialized");

    Density &rho = this->density;
    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &X = *this->orbitals_x;
    OrbitalVector &Y = *this->orbitals_y;

    Timer timer;
    density::compute(prec, rho, Phi, X, Y, DensityType::Total);
    print_utils::qmfunction(3, "Compute global density", rho, timer);
}

void CoulombPotentialD2::setupLocalDensity(double prec) {
    if (hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_x == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_y == nullptr) MSG_ERROR("Orbitals not initialized");

    Density &rho = this->density;
    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &X = *this->orbitals_x;
    OrbitalVector &Y = *this->orbitals_y;

    Timer timer;
    density::compute_local(prec, rho, Phi, X, Y, DensityType::Total);
    print_utils::qmfunction(3, "Compute local density", rho, timer);
}

} // namespace mrchem
