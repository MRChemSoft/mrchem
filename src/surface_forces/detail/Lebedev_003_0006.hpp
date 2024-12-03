
#pragma once

#include <array>
#include <tuple>

#include <Eigen/Core>

#include "lebedev_utils.hpp"
#include "meta_utils.hpp"

namespace detail {
class Lebedev_003_0006 final
{
private:
  /** Number of points in A1 block */
  static constexpr auto N_a1 = 6;
  /** Quadrature weight in A1 block */
  static constexpr auto a1 = 0.16666666666666666;

public:
  using Weights = Eigen::VectorXd;
  using Points = Eigen::Matrix3Xd;

  static auto quadrature() -> std::tuple<Weights, Points>
  {
    Weights ws = Weights::Zero(6);
    Points xs = Points::Zero(3, 6);

    {
      auto ws_ = ws(Eigen::seqN(Eigen::fix<0>, Eigen::fix<N_a1>));
      auto xs_ = xs(Eigen::all, Eigen::seqN(Eigen::fix<0>, Eigen::fix<N_a1>));
      std::tie(ws_, xs_) = detail::a1(a1);
    }

    return { ws, xs };
  }
};
} // namespace detail
