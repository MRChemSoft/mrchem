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

#include <nlohmann/json.hpp>

namespace mrchem {
namespace mrenv {

nlohmann::json fetch_json(int argc, char **argv);
void dump_json(const nlohmann::json &json_inp, const nlohmann::json &json_out);

void initialize(const nlohmann::json &json_inp);
void finalize(double wt);

} // namespace mrenv
namespace detail {
std::string remove_extension(const std::string &fname);
bool all_success(const nlohmann::json &json_out);
} // namespace detail
} // namespace mrchem
