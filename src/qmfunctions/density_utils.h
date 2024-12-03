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

#include "mrchem.h"
#include "Orbital.h"
#include "Density.h"

namespace mrchem {
namespace density {

void allreduce_density(double prec, Density &rho_tot, Density &rho_loc);
void compute(double prec, Density &rho, mrcpp::GaussExp<3> &dens_exp);
void compute(double prec, Density &rho, OrbitalVector &Phi, DensityType spin);
void compute(double prec, Density &rho, OrbitalVector &Phi, OrbitalVector &X, OrbitalVector &Y, DensityType spin);
void compute_local(double prec, Density &rho, OrbitalVector &Phi, DensityType spin);
void compute_local(double prec, Density &rho, OrbitalVector &Phi, OrbitalVector &X, OrbitalVector &Y, DensityType spin);

/**
 * @brief Reads atomic density data from a file
 * @param path The path to the file containing the data
 * @param rGrid The grid for the radial distances
 * @param rhoGrid The grid for the atomic density values
 */
void readAtomicDensity(const std::string path, Eigen::VectorXd &rGrid, Eigen::VectorXd &rhoGrid);

} // namespace density
} // namespace mrchem
