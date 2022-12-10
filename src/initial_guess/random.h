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

#include <Eigen/Dense>

#include "qmfunctions/qmfunction_fwd.h"

/** @file random.h
 *
 * @brief Module for generating response calculations initial guess from
 * application of a linear combination, with random coefficients, of products of
 * monomials in x, y, z.
 */

namespace mrchem {
namespace initial_guess {
namespace random {

bool setup(OrbitalVector &XorY, OrbitalVector &Phi, double prec, int seed);

bool apply_random(OrbitalVector &XorY, const Eigen::Vector3d &coeffs, double prec);

} // namespace random
} // namespace initial_guess
} // namespace mrchem
