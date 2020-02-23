/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "qmfunctions/qmfunction_fwd.h"

/** @file sad.h
 *
 * @brief Module for generating initial guess as superposition of atomic densities
 *
 * The initial_guess::sad namespace provides functionality to setup an initial
 * guess from hydrogen functions and a superposition of atomic densities.
 */

namespace mrchem {
class Nuclei;

namespace initial_guess {
namespace sad {

bool setup(OrbitalVector &Phi, double prec, const Nuclei &nucs, int zeta);

} // namespace sad
} // namespace initial_guess
} // namespace mrchem
