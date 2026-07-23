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
#include "qmfunctions/Orbital.h"
#include "qmfunctions/qmfunction_fwd.h"
#include <unsupported/Eigen/CXX11/Tensor>

namespace mrchem {
/** @brief precomputes all the <ik|1/r12|jl> 2-electron integrals
 *
 *  @param[in] Phi input orbitals
 *
 * The ij/r12 potentials are first (pre)computed for all the orbitals j<i, then multiplied by kl
 */
Eigen::Tensor<std::complex<double>, 4> calc_2elintegrals(double prec, OrbitalVector &Phi);
/** @brief computes Int(phi_i^dag*phi_j/|r-r'|)
 *
 *  \param[in] phi_i orbital to be conjugated and multiplied by phi_j
 *  \param[in] phi_j orbital to be multiplied by phi_i^dag
 *  \param[in] rho is normally the density. It is used to screen the final output
 *  \param[out] V_ij is the resulting potential function
  *
 * Computes the product of complex conjugate of phi_i and phi_j,
 * then applies the Poisson operator. The result is given in V_ij.
 */
void ij_r12(double prec, Orbital rho, Orbital phi_i, Orbital phi_j, Orbital &V_ij) ;
}
