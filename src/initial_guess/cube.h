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

#include <highfive/H5Easy.hpp>
#include <nlohmann/json.hpp>

#include "analyticfunctions/CUBEfunction.h"
#include "qmfunctions/qmfunction_fwd.h"

/** @file cube.h
 *
 * @brief Module for generating initial guess from cube file
 *
 */

namespace mrchem {
namespace initial_guess {
namespace cube {
bool setup(OrbitalVector &Phi, double prec, const std::string &file_p, const std::string &file_a, const std::string &file_b);
bool project_mo(OrbitalVector &Phi, double prec, const std::string &mo_file);
std::vector<mrchem::CUBEfunction> getCUBEFunction(const nlohmann::json &json_cube);
} // namespace cube

namespace hdf5 {
bool setup(OrbitalVector &Phi, double prec, const std::string &h5fname);
bool project_mo(OrbitalVector &Phi, double prec, const H5Easy::File &h5file, const std::string &spin_lbl);
std::vector<mrchem::CUBEfunction> getCUBEFunction(const H5Easy::File &h5file, const std::string &spin_lbl);
} // namespace hdf5
} // namespace initial_guess
} // namespace mrchem
