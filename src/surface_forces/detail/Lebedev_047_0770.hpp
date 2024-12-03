
#pragma once

#include <array>
#include <tuple>

#include <Eigen/Core>

#include "lebedev_utils.hpp"
#include "meta_utils.hpp"

namespace detail {
class Lebedev_047_0770 final
{
private:
  /** Number of points in A1 block */
  static constexpr auto N_a1 = 6;
  /** Quadrature weight in A1 block */
  static constexpr auto a1 = 0.0002192942088182175;

  /** Number of points in A2 block */
  static constexpr auto N_a2 = 12;
  /** Quadrature weight in A2 block */
  static constexpr auto a2 = 0.001436433617319091;

  /** Number of points in A3 block */
  static constexpr auto N_a3 = 8;
  /** Quadrature weight in A3 block */
  static constexpr auto a3 = 0.001421940344335874;

  /** Number of points in each PQ0 block */
  static constexpr auto N_pq0 = 24;
  /** Quadrature weight(s) and point(s) in PQ0 block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as \f$(p, q,
   * 0) = (\sin(x * \pi), \cos(x * \pi), 0)\f$
   */
  static constexpr std::array<std::array<double, 2>, 3> pq0 = { { { { 0.0009254401499865399, 0.04621738571378741 } },
                                                                  { { 0.001250239995053508, 0.1100975833825021 } },
                                                                  { { 0.001394365843329228, 0.1791538510007847 } } } };

  /** Number of points in each LLM block */
  static constexpr auto N_llm = 24;
  /** Quadrature weight(s) and point(s) in LLM block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(l, l, m) = (\frac{\sin(x * \pi)}{\sqrt{2}}, \frac{\sin(x * \pi)}{\sqrt{2}}, \cos(x * \pi))\f$
   */
  static constexpr std::array<std::array<double, 2>, 10> llm = { { { { 0.0006798123511050234, 0.02292026702769342 } },
                                                                   { { 0.0009913184235294894, 0.0555702172024949 } },
                                                                   { { 0.001180207833238949, 0.09254023339464029 } },
                                                                   { { 0.001296599602080922, 0.1319393297241756 } },
                                                                   { { 0.001365871427428312, 0.1730026707517375 } },
                                                                   { { 0.001402988604775323, 0.2154115731074323 } },
                                                                   { { 0.001418645563595606, 0.2590829375936784 } },
                                                                   { { 0.001421376741851662, 0.3505993811839581 } },
                                                                   { { 0.001423996475490966, 0.3988264338132601 } },
                                                                   { { 0.001431554042178568, 0.4488003348564934 } } } };

  /** Number of points in each RSW block */
  static constexpr auto N_rsw = 48;
  /** Quadrature weight(s) and point(s) in RSW block.
   *
   * Each element is a {weight, x, y} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(r, s, w) = (\sin(x * \pi)\cos(y * \pi), \sin(x * \pi)\sin(y * \pi), \cos(x*\pi))\f$
   */
  static constexpr std::array<std::array<double, 3>, 9> rsw = {
    { { { 0.00112708909467175, 0.07896604124950307, 0.09126444167498717 } },
      { { 0.001345753760910667, 0.1553094137600919, 0.1608217314832364 } },
      { { 0.001424957283316787, 0.215550404070437, 0.0408824237566492 } },
      { { 0.001261523341237749, 0.115923272601319, 0.1352197040723114 } },
      { { 0.001392547106052695, 0.1836874558196702, 0.09300622918241049 } },
      { { 0.001418761677877655, 0.2238040826792707, 0.119230464563211 } },
      { { 0.001338366684479552, 0.1455361124992715, 0.05591508238246502 } },
      { { 0.00139370086267613, 0.1964019716288602, 0.1775706977927038 } },
      { { 0.00141591475746693, 0.2389135816551385, 0.1893661429124828 } } }
  };

public:
  using Weights = Eigen::VectorXd;
  using Points = Eigen::Matrix3Xd;

  static auto quadrature() -> std::tuple<Weights, Points>
  {
    Weights ws = Weights::Zero(770);
    Points xs = Points::Zero(3, 770);

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
