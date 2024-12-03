
#pragma once

#include <array>
#include <tuple>

#include <Eigen/Core>

#include "lebedev_utils.hpp"
#include "meta_utils.hpp"

namespace detail {
class Lebedev_035_0434 final
{
private:
  /** Number of points in A1 block */
  static constexpr auto N_a1 = 6;
  /** Quadrature weight in A1 block */
  static constexpr auto a1 = 0.0005265897968228898;

  /** Number of points in A2 block */
  static constexpr auto N_a2 = 12;
  /** Quadrature weight in A2 block */
  static constexpr auto a2 = 0.002548219972002607;

  /** Number of points in A3 block */
  static constexpr auto N_a3 = 8;
  /** Quadrature weight in A3 block */
  static constexpr auto a3 = 0.002512317418927299;

  /** Number of points in each PQ0 block */
  static constexpr auto N_pq0 = 24;
  /** Quadrature weight(s) and point(s) in PQ0 block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as \f$(p, q,
   * 0) = (\sin(x * \pi), \cos(x * \pi), 0)\f$
   */
  static constexpr std::array<std::array<double, 2>, 2> pq0 = { { { { 0.002417442375638987, 0.1563228953888757 } },
                                                                  { { 0.001910951282179527, 0.06743512918521113 } } } };

  /** Number of points in each LLM block */
  static constexpr auto N_llm = 24;
  /** Quadrature weight(s) and point(s) in LLM block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(l, l, m) = (\frac{\sin(x * \pi)}{\sqrt{2}}, \frac{\sin(x * \pi)}{\sqrt{2}}, \cos(x * \pi))\f$
   */
  static constexpr std::array<std::array<double, 2>, 7> llm = { { { { 0.002530403801186355, 0.4317916211421352 } },
                                                                  { { 0.002014279020918488, 0.08075932475954789 } },
                                                                  { { 0.002501725168402938, 0.2445925526145315 } },
                                                                  { { 0.002513267174597568, 0.3663265069319693 } },
                                                                  { { 0.002302694782227404, 0.1326057309579972 } },
                                                                  { { 0.00146249562159454, 0.0341337299597072 } },
                                                                  { { 0.002445373437312978, 0.1874342228194836 } } } };

  /** Number of points in each RSW block */
  static constexpr auto N_rsw = 48;
  /** Quadrature weight(s) and point(s) in RSW block.
   *
   * Each element is a {weight, x, y} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(r, s, w) = (\sin(x * \pi)\cos(y * \pi), \sin(x * \pi)\sin(y * \pi), \cos(x*\pi))\f$
   */
  static constexpr std::array<std::array<double, 3>, 4> rsw = {
    { { { 0.002416930044324782, 0.1647978238635747, 0.1362863107481372 } },
      { { 0.002512236854563495, 0.2048710360686823, 0.05695446177756676 } },
      { { 0.002496644054553092, 0.2193835459268291, 0.1623274676288571 } },
      { { 0.002236607760437855, 0.1134255321813034, 0.09180065418907327 } } }
  };

public:
  using Weights = Eigen::VectorXd;
  using Points = Eigen::Matrix3Xd;

  static auto quadrature() -> std::tuple<Weights, Points>
  {
    Weights ws = Weights::Zero(434);
    Points xs = Points::Zero(3, 434);

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

    // number of pq0 values
    constexpr auto pq0_sz = pq0.size();
    detail::constexpr_for<0, pq0_sz, 1>([&xs, &ws](auto i) {
      auto ws_ = ws(Eigen::seqN(Eigen::fix<+N_a1 + N_a2 + N_a3 + i * N_pq0>, Eigen::fix<N_pq0>));
      auto xs_ = xs(Eigen::all, Eigen::seqN(Eigen::fix<+N_a1 + N_a2 + N_a3 + i * N_pq0>, Eigen::fix<N_pq0>));
      std::tie(ws_, xs_) = detail::pq0(pq0[i][0], pq0[i][1]);
    });

    // number of llm values
    constexpr auto llm_sz = llm.size();
    detail::constexpr_for<0, llm_sz, 1>([&xs, &ws](auto i) {
      auto ws_ = ws(Eigen::seqN(Eigen::fix<+N_a1 + N_a2 + N_a3 + pq0.size() * N_pq0 + i * N_llm>, Eigen::fix<N_llm>));
      auto xs_ = xs(Eigen::all,
                    Eigen::seqN(Eigen::fix<+N_a1 + N_a2 + N_a3 + pq0.size() * N_pq0 + i * N_llm>, Eigen::fix<N_llm>));
      std::tie(ws_, xs_) = detail::llm(llm[i][0], llm[i][1]);
    });

    // number of rsw values
    constexpr auto rsw_sz = rsw.size();
    detail::constexpr_for<0, rsw_sz, 1>([&xs, &ws](auto i) {
      auto ws_ = ws(Eigen::seqN(Eigen::fix<+N_a1 + N_a2 + N_a3 + pq0.size() * N_pq0 + llm.size() * N_llm + i * N_rsw>,
                                Eigen::fix<N_rsw>));
      auto xs_ = xs(Eigen::all,
                    Eigen::seqN(Eigen::fix<+N_a1 + N_a2 + N_a3 + pq0.size() * N_pq0 + llm.size() * N_llm + i * N_rsw>,
                                Eigen::fix<N_rsw>));
      std::tie(ws_, xs_) = detail::rsw(rsw[i][0], rsw[i][1], rsw[i][2]);
    });

    return { ws, xs };
  }
};
} // namespace detail
