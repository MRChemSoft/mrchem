
#pragma once

#include <array>
#include <tuple>

#include <Eigen/Core>

#include "lebedev_utils.hpp"
#include "meta_utils.hpp"

namespace detail {
class Lebedev_053_0974 final
{
private:
  /** Number of points in A1 block */
  static constexpr auto N_a1 = 6;
  /** Quadrature weight in A1 block */
  static constexpr auto a1 = 0.000143829419053;

  /** Number of points in A3 block */
  static constexpr auto N_a3 = 8;
  /** Quadrature weight in A3 block */
  static constexpr auto a3 = 0.001125772288287;

  /** Number of points in each PQ0 block */
  static constexpr auto N_pq0 = 24;
  /** Quadrature weight(s) and point(s) in PQ0 block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as \f$(p, q,
   * 0) = (\sin(x * \pi), \cos(x * \pi), 0)\f$
   */
  static constexpr std::array<std::array<double, 2>, 4> pq0 = { { { { 0.000682336792711, 0.03946600009801534 } },
                                                                  { { 0.000945415816045, 0.09501245088959068 } },
                                                                  { { 0.001074429975386, 0.15566509596376388 } },
                                                                  { { 0.001129300086569, 0.2183643706356081 } } } };

  /** Number of points in each LLM block */
  static constexpr auto N_llm = 24;
  /** Quadrature weight(s) and point(s) in LLM block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(l, l, m) = (\frac{\sin(x * \pi)}{\sqrt{2}}, \frac{\sin(x * \pi)}{\sqrt{2}}, \cos(x * \pi))\f$
   */
  static constexpr std::array<std::array<double, 2>, 12> llm = { { { { 0.000494802934195, 0.019337017127159052 } },
                                                                   { { 0.000735799010913, 0.04750701000375154 } },
                                                                   { { 0.00088891327713, 0.07960603068581106 } },
                                                                   { { 0.000988834783892, 0.1139520045627366 } },
                                                                   { { 0.001053299681709, 0.1498166059474992 } },
                                                                   { { 0.001092778807015, 0.18685324600959904 } },
                                                                   { { 0.001114389394063, 0.22491326381309878 } },
                                                                   { { 0.001123724788052, 0.26397058054454303 } },
                                                                   { { 0.001125239325244, 0.34538559044616385 } },
                                                                   { { 0.001126153271816, 0.38800992457940353 } },
                                                                   { { 0.001130286931124, 0.4320172490827723 } },
                                                                   { { 0.001134986534364, 0.4771817066440508 } } } };

  /** Number of points in each RSW block */
  static constexpr auto N_rsw = 48;
  /** Quadrature weight(s) and point(s) in RSW block.
   *
   * Each element is a {weight, x, y} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(r, s, w) = (\sin(x * \pi)\cos(y * \pi), \sin(x * \pi)\sin(y * \pi), \cos(x*\pi))\f$
   */
  static constexpr std::array<std::array<double, 3>, 12> rsw = {
    { { { 0.00084368845009, 0.0678427744589942, 0.0911399534108118 } },
      { { 0.001075255720449, 0.15950256677371893, 0.09245458231700068 } },
      { { 0.001108577236864, 0.20750107737350024, 0.18879250210965065 } },
      { { 0.000956647532378, 0.10003743318687472, 0.13495092827433242 } },
      { { 0.001080663250717, 0.17037391906706956, 0.17706122532834084 } },
      { { 0.001126797131196, 0.22169791057730812, 0.07103466827059002 } },
      { { 0.001022568715358, 0.12604674065566943, 0.0556505520034843 } },
      { { 0.001108960267713, 0.1877588469112349, 0.04049681224647043 } },
      { { 0.001122790653436, 0.2311961883820272, 0.13768263611849016 } },
      { { 0.001032401847117, 0.13444579383088165, 0.16042026040349594 } },
      { { 0.001107249382284, 0.19466826141718838, 0.118434013056699 } },
      { { 0.00112178004852, 0.24570106099668815, 0.19748919518870767 } } }
  };

public:
  using Weights = Eigen::VectorXd;
  using Points = Eigen::Matrix3Xd;

  static auto quadrature() -> std::tuple<Weights, Points>
  {
    Weights ws = Weights::Zero(974);
    Points xs = Points::Zero(3, 974);

    {
      auto ws_ = ws(Eigen::seqN(Eigen::fix<0>, Eigen::fix<N_a1>));
      auto xs_ = xs(Eigen::all, Eigen::seqN(Eigen::fix<0>, Eigen::fix<N_a1>));
      std::tie(ws_, xs_) = detail::a1(a1);
    }

    {
      auto ws_ = ws(Eigen::seqN(Eigen::fix<+N_a1>, Eigen::fix<N_a3>));
      auto xs_ = xs(Eigen::all, Eigen::seqN(Eigen::fix<+N_a1>, Eigen::fix<N_a3>));
      std::tie(ws_, xs_) = detail::a3(a3);
    }

    // number of pq0 values
    constexpr auto pq0_sz = pq0.size();
    detail::constexpr_for<0, pq0_sz, 1>([&xs, &ws](auto i) {
      auto ws_ = ws(Eigen::seqN(Eigen::fix<+N_a1 + N_a3 + i * N_pq0>, Eigen::fix<N_pq0>));
      auto xs_ = xs(Eigen::all, Eigen::seqN(Eigen::fix<+N_a1 + N_a3 + i * N_pq0>, Eigen::fix<N_pq0>));
      std::tie(ws_, xs_) = detail::pq0(pq0[i][0], pq0[i][1]);
    });

    // number of llm values
    constexpr auto llm_sz = llm.size();
    detail::constexpr_for<0, llm_sz, 1>([&xs, &ws](auto i) {
      auto ws_ = ws(Eigen::seqN(Eigen::fix<+N_a1 + N_a3 + pq0.size() * N_pq0 + i * N_llm>, Eigen::fix<N_llm>));
      auto xs_ =
        xs(Eigen::all, Eigen::seqN(Eigen::fix<+N_a1 + N_a3 + pq0.size() * N_pq0 + i * N_llm>, Eigen::fix<N_llm>));
      std::tie(ws_, xs_) = detail::llm(llm[i][0], llm[i][1]);
    });

    // number of rsw values
    constexpr auto rsw_sz = rsw.size();
    detail::constexpr_for<0, rsw_sz, 1>([&xs, &ws](auto i) {
      auto ws_ = ws(
        Eigen::seqN(Eigen::fix<+N_a1 + N_a3 + pq0.size() * N_pq0 + llm.size() * N_llm + i * N_rsw>, Eigen::fix<N_rsw>));
      auto xs_ = xs(
        Eigen::all,
        Eigen::seqN(Eigen::fix<+N_a1 + N_a3 + pq0.size() * N_pq0 + llm.size() * N_llm + i * N_rsw>, Eigen::fix<N_rsw>));
      std::tie(ws_, xs_) = detail::rsw(rsw[i][0], rsw[i][1], rsw[i][2]);
    });

    return { ws, xs };
  }
};
} // namespace detail
