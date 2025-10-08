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
#include "qmfunctions/qmfunction_fwd.h"
#include "tensor/tensor_fwd.h"

#include "qmoperators/one_electron/KineticOperator.h"

/** @file core.h
 *
 * @brief Module for generating initial guess of hydrogen functions
 *
 * The initial_guess::core namespace provides functionality to setup an
 * initial guess of hydrogen eigenfunctions.
 */

namespace mrchem {
class Nuclei;

namespace initial_guess {
namespace core {

bool setup(OrbitalVector &Phi, double prec, const Nuclei &nucs, int zeta);
void project_ao(OrbitalVector &Phi, double prec, const Nuclei &nucs, int zeta);
void rotate_orbitals(OrbitalVector &Psi, double prec, ComplexMatrix &U, OrbitalVector &Phi);
ComplexMatrix diagonalize(OrbitalVector &Phi, MomentumOperator &T, RankZeroOperator &V);
ComplexMatrix diagonalize(OrbitalVector &Phi, ComplexMatrix &t_tilde, RankZeroOperator &V);

} // namespace core
} // namespace initial_guess
} // namespace mrchem
