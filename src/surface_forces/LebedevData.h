#pragma once

#include <tuple>

#include <Eigen/Core>

/** Lebedev-Laikov partition with at most \f$N\f$ points.
 *
 * @param[in] N
 * @return the weights and points.
 */
auto
lebedev(size_t N) -> std::tuple<Eigen::VectorXd, Eigen::Matrix3Xd>;