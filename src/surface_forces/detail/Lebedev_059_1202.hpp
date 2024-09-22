
#pragma once

#include <array>
#include <tuple>

#include <Eigen/Core>

#include "lebedev_utils.hpp"
#include "meta_utils.hpp"

namespace detail {
class Lebedev_059_1202 final
{
private:
  /** Number of points in A1 block */
  static constexpr auto N_a1 = 6;
  /** Quadrature weight in A1 block */
  static constexpr auto a1 = 0.000110518923327;

  /** Number of points in A2 block */
  static constexpr auto N_a2 = 12;
  /** Quadrature weight in A2 block */
  static constexpr auto a2 = 0.000920523273809;

  /** Number of points in A3 block */
  static constexpr auto N_a3 = 8;
  /** Quadrature weight in A3 block */
  static constexpr auto a3 = 0.000913315978644;

  /** Number of points in each PQ0 block */
  static constexpr auto N_pq0 = 24;
  /** Quadrature weight(s) and point(s) in PQ0 block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as \f$(p, q,
   * 0) = (\sin(x * \pi), \cos(x * \pi), 0)\f$
   */
  static constexpr std::array<std::array<double, 2>, 4> pq0 = { { { { 0.000517697731297, 0.03420075260565257 } },
                                                                  { { 0.00073311436821, 0.08313162364341634 } },
                                                                  { { 0.000846323283638, 0.13701491459108198 } },
                                                                  { { 0.000903112269425, 0.19307102433703996 } } } };

  /** Number of points in each LLM block */
  static constexpr auto N_llm = 24;
  /** Quadrature weight(s) and point(s) in LLM block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(l, l, m) = (\frac{\sin(x * \pi)}{\sqrt{2}}, \frac{\sin(x * \pi)}{\sqrt{2}}, \cos(x * \pi))\f$
   */
  static constexpr std::array<std::array<double, 2>, 13> llm = { { { { 0.000369042189802, 0.016720424160189817 } },
                                                                   { { 0.000560399092868, 0.04126017303586406 } },
                                                                   { { 0.000686529762928, 0.0694729706214065 } },
                                                                   { { 0.000772033855115, 0.09980353185791829 } },
                                                                   { { 0.000830154595889, 0.13155235170779025 } },
                                                                   { { 0.000868669255018, 0.16436581638357306 } },
                                                                   { { 0.000892707628585, 0.1980657851787268 } },
                                                                   { { 0.000906082023857, 0.23257851288553447 } },
                                                                   { { 0.000911977725494, 0.26790167541606985 } },
                                                                   { { 0.00091287201386, 0.3412238232237018 } },
                                                                   { { 0.000913071493569, 0.3794156291080754 } },
                                                                   { { 0.000915287378455, 0.4187196099185521 } },
                                                                   { { 0.000918743627432, 0.45903955700964827 } } } };

  /** Number of points in each RSW block */
  static constexpr auto N_rsw = 48;
  /** Quadrature weight(s) and point(s) in RSW block.
   *
   * Each element is a {weight, x, y} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(r, s, w) = (\sin(x * \pi)\cos(y * \pi), \sin(x * \pi)\sin(y * \pi), \cos(x*\pi))\f$
   */
  static constexpr std::array<std::array<double, 3>, 16> rsw = {
    { { { 0.000648577845316, 0.059125048710861904, 0.09104538585311528 } },
      { { 0.000743503091098, 0.08753554058354115, 0.13476625655108154 } },
      { { 0.000799852789184, 0.11063719205345883, 0.05547084389337896 } },
      { { 0.000810173149747, 0.117990254123381, 0.16014183229915946 } },
      { { 0.000848338957459, 0.14034322517430972, 0.09207018989157909 } },
      { { 0.000855629925731, 0.1498356879851596, 0.17669748867405466 } },
      { { 0.000880320867974, 0.1656336092513369, 0.040224891242256824 } },
      { { 0.000881104818243, 0.1715841534183696, 0.11785881091491886 } },
      { { 0.000885028234127, 0.18274525247006784, 0.18836275020540447 } },
      { { 0.000902134229904, 0.19589548091013856, 0.07047600592771162 } },
      { { 0.000901009167711, 0.20401067271108056, 0.1369616236187671 } },
      { { 0.000902269293843, 0.2165602520928162, 0.19702704234601376 } },
      { { 0.000915801617469, 0.22216394131587347, 0.03187470178277977 } },
      { { 0.000913157800319, 0.22746279397069277, 0.09394027570063196 } },
      { { 0.000910781357948, 0.23744932011918268, 0.15164004637018952 } },
      { { 0.000910576025897, 0.25122346068669016, 0.20369796009147373 } } }
  };

public:
  using Weights = Eigen::VectorXd;
  using Points = Eigen::Matrix3Xd;

  static auto quadrature() -> std::tuple<Weights, Points>
  {
    Weights ws = Weights::Zero(1202);
    Points xs = Points::Zero(3, 1202);

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
