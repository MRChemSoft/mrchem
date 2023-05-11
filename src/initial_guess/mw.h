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

#include <string>

#include "qmfunctions/qmfunction_fwd.h"

/** @file mw.h
 *
 * @brief Module for generating initial guess from previous MW calculations
 *
 * The initial_guess::mw namespace provides functionality to setup an initial
 * guess from orbitals stored in previous MW calculations.
 */

namespace mrchem {
namespace initial_guess {
namespace mw {

bool setup(OrbitalVector &Phi, double prec, const std::string &file_p, const std::string &file_a, const std::string &file_b);

} // namespace mw
} // namespace initial_guess
} // namespace mrchem
