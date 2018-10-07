/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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

/** Notes on vectors:
 * The OrbitalVector (std::vector<Orbital>) is conceptually different from the
 * QMFunctionVector (std::vector<QMFunction>). The former should be a collection
 * of orbitals defining a proper Slater determinant, while the latter is any
 * arbitrary collection of functions (could also be orbitals) used e.g. when
 * adding several orbitals into one. The difference is important when MPI is
 * considered, as OrbitalVectors are distributed while QMFunctionVectors are
 * not. This is reflected in the functionality that is available for the
 * different vectors, where OrbitalVectors should always be treated as a whole,
 * e.g. in an orbital rotation Psi = orbital::rotate(U, Phi), where U is a
 * unitary matrix. This operation is MPI safe and all necessary communication
 * happen automatically under the hood. This is opposed to the similar linear
 * combination qmfunction::linear_combination(func_out, coefs_vec, func_vec)
 * where no MPI is considered.
 */

namespace mrchem {

class QMFunction;
using QMFunctionVector = std::vector<QMFunction>;

class Orbital;
using OrbitalChunk = std::vector<std::tuple<int, Orbital>>;
using OrbitalVector = std::vector<Orbital>;

class Density;

} //namespace mrchem
