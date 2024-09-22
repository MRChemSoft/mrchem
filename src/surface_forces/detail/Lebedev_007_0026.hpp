
#pragma once

#include <array>
#include <tuple>

#include <Eigen/Core>

#include "lebedev_utils.hpp"
#include "meta_utils.hpp"

namespace detail {
class Lebedev_007_0026 final
{
private:
  /** Number of points in A1 block */
  static constexpr auto N_a1 = 6;
  /** Quadrature weight in A1 block */
  static constexpr auto a1 = 0.04761904761904772;

  /** Number of points in A2 block */
  static constexpr auto N_a2 = 12;
  /** Quadrature weight in A2 block */
  static constexpr auto a2 = 0.03809523809523815;

  /** Number of points in A3 block */
  static constexpr auto N_a3 = 8;
  /** Quadrature weight in A3 block */
  static constexpr auto a3 = 0.03214285714285716;

public:
  using Weights = Eigen::VectorXd;
  using Points = Eigen::Matrix3Xd;

  static auto quadrature() -> std::tuple<Weights, Points>
  {
    Weights ws = Weights::Zero(26);
    Points xs = Points::Zero(3, 26);

    {
      auto ws_ = ws(Eigen::seqN(Eigen::fix<0>, Eigen::fix<N_a1>));
      auto xs_ = xs(Eigen::all, Eigen::seqN(Eigen::fix<0>, Eigen::fix<N_a1>));
      std::tie(ws_, xs_) = detail::a1(a1);
    }

    {
      auto ws_ = ws(Eigen::seqN(Eigen::fix<+N_a1>, Eigen::fix<N_a2>));
      auto xs_ = xs(Eigen::all, Eigen::seqN(Eigen::fix<+N_a1>, Eigen::fix<N_a2>));
      std::tie(ws_, xs_) = detail::a2(a2);
    }

    {
      auto ws_ = ws(Eigen::seqN(Eigen::fix<+N_a1 + N_a2>, Eigen::fix<N_a3>));
      auto xs_ = xs(Eigen::all, Eigen::seqN(Eigen::fix<+N_a1 + N_a2>, Eigen::fix<N_a3>));
      std::tie(ws_, xs_) = detail::a3(a3);
    }

    return { ws, xs };
  }
};
} // namespace detail
