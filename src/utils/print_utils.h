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

#include <nlohmann/json.hpp>

namespace mrchem {
namespace print_utils {

void centered_text(int level, const std::string &txt);
void headline(int level, const std::string &txt);
void text(int level, const std::string &txt, const std::string &val, bool ralign = false);
void json(int level, const nlohmann::json &j, bool ralign = false);
void coord(int level, const std::string &txt, const mrcpp::Coord<3> &val, int p = -1, bool s = false);
void scalar(int level, const std::string &txt, double val, const std::string &unit = "", int p = -1, bool s = false);
void vector(int level, const std::string &txt, const DoubleVector &val, int p = -1, bool s = false);
void matrix(int level, const std::string &txt, const DoubleMatrix &val, int p = -1, bool s = false);
void qmfunction(int level, const std::string &txt, const mrcpp::ComplexFunction &func, mrcpp::Timer &t);
std::string dbl_to_str(double d, int p, bool sci);
std::vector<double> eigen_to_vector(const DoubleVector &inp, double thrs);
std::vector<double> eigen_to_vector(const DoubleMatrix &inpm, double thrs);
void ltrim(std::string &s);
void rtrim(std::string &s);
void trim(std::string &s);
std::string ltrim_copy(std::string s);
std::string rtrim_copy(std::string s);
std::string trim_copy(std::string s);

} // namespace print_utils
} // namespace mrchem
