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

#include <tuple>

#include <Eigen/Core>


namespace detail {
auto
a1(double w) -> std::tuple<Eigen::Matrix<double, 6, 1>, Eigen::Matrix<double, 3, 6>>;

auto
a2(double w) -> std::tuple<Eigen::Matrix<double, 12, 1>, Eigen::Matrix<double, 3, 12>>;

auto
a3(double w) -> std::tuple<Eigen::Matrix<double, 8, 1>, Eigen::Matrix<double, 3, 8>>;

auto
pq0(double w, double x) -> std::tuple<Eigen::Matrix<double, 24, 1>, Eigen::Matrix<double, 3, 24>>;

auto
llm(double w, double x) -> std::tuple<Eigen::Matrix<double, 24, 1>, Eigen::Matrix<double, 3, 24>>;

auto
rsw(double w, double x, double y) -> std::tuple<Eigen::Matrix<double, 48, 1>, Eigen::Matrix<double, 3, 48>>;

auto
n_points_to_degree(size_t N) -> size_t;
} // namespace detail

