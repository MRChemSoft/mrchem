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

/** @file math_utils.h
 *
 * @brief Collection of stand-alone math related functions
 *
 */

namespace mrchem {
namespace math_utils {

double calc_distance(const mrcpp::Coord<3> &a, const mrcpp::Coord<3> &b);

DoubleVector init_nan(int I);
DoubleMatrix init_nan(int I, int J);

DoubleMatrix read_matrix_file(const std::string &file);
DoubleMatrix skew_matrix_exp(const DoubleMatrix &A);
ComplexMatrix hermitian_matrix_pow(const ComplexMatrix &A, double b);
ComplexMatrix diagonalize_hermitian_matrix(const ComplexMatrix &A, DoubleVector &diag);
void diagonalize_block(ComplexMatrix &M, ComplexMatrix &U, int nstart, int nsize);
double logsumexp(const Eigen::VectorXd &x);

} // namespace math_utils
} // namespace mrchem
