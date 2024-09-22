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

