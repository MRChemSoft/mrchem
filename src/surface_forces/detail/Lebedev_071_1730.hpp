
#pragma once

#include <array>
#include <tuple>

#include <Eigen/Core>

#include "lebedev_utils.hpp"
#include "meta_utils.hpp"

namespace detail {
class Lebedev_071_1730 final
{
private:
  /** Number of points in A1 block */
  static constexpr auto N_a1 = 6;
  /** Quadrature weight in A1 block */
  static constexpr auto a1 = 6.3090494374e-05;

  /** Number of points in A2 block */
  static constexpr auto N_a2 = 12;
  /** Quadrature weight in A2 block */
  static constexpr auto a2 = 0.000639828770557;

  /** Number of points in A3 block */
  static constexpr auto N_a3 = 8;
  /** Quadrature weight in A3 block */
  static constexpr auto a3 = 0.000635718507353;

  /** Number of points in each PQ0 block */
  static constexpr auto N_pq0 = 24;
  /** Quadrature weight(s) and point(s) in PQ0 block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as \f$(p, q,
   * 0) = (\sin(x * \pi), \cos(x * \pi), 0)\f$
   */
  static constexpr std::array<std::array<double, 2>, 5> pq0 = { { { { 0.000318691344995, 0.02662567388107946 } },
                                                                  { { 0.000467802855859, 0.06572409006333814 } },
                                                                  { { 0.00055388296976, 0.1094192803396822 } },
                                                                  { { 0.000604447590719, 0.15535139330769052 } },
                                                                  { { 0.000631357510351, 0.20242837200985242 } } } };

  /** Number of points in each LLM block */
  static constexpr auto N_llm = 24;
  /** Quadrature weight(s) and point(s) in LLM block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(l, l, m) = (\frac{\sin(x * \pi)}{\sqrt{2}}, \frac{\sin(x * \pi)}{\sqrt{2}}, \cos(x * \pi))\f$
   */
  static constexpr std::array<std::array<double, 2>, 16> llm = { { { { 0.000222120716219, 0.0128821951216131 } },
                                                                   { { 0.000347578402229, 0.032207731604644504 } },
                                                                   { { 0.000435074244359, 0.054701957361929614 } },
                                                                   { { 0.000497856913652, 0.07907872541215494 } },
                                                                   { { 0.0005435036222, 0.10472423247884396 } },
                                                                   { { 0.000576591338822, 0.1313033365858972 } },
                                                                   { { 0.000600120035923, 0.15862323583036297 } },
                                                                   { { 0.000616217817272, 0.1865740566052706 } },
                                                                   { { 0.000626521815244, 0.21509963472632135 } },
                                                                   { { 0.000632398716097, 0.2441820363305423 } },
                                                                   { { 0.000635076785154, 0.27383276117335115 } },
                                                                   { { 0.00063543627753, 0.33499578715110123 } },
                                                                   { { 0.000635230246271, 0.3666178126756868 } },
                                                                   { { 0.000635811788142, 0.39899550033190523 } },
                                                                   { { 0.000637310159031, 0.43211887236649743 } },
                                                                   { { 0.000639042896137, 0.46587318374717057 } } } };

  /** Number of points in each RSW block */
  static constexpr auto N_rsw = 48;
  /** Quadrature weight(s) and point(s) in RSW block.
   *
   * Each element is a {weight, x, y} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(r, s, w) = (\sin(x * \pi)\cos(y * \pi), \sin(x * \pi)\sin(y * \pi), \cos(x*\pi))\f$
   */
  static constexpr std::array<std::array<double, 3>, 25> rsw = {
    { { { 0.000407862643186, 0.04644429609109567, 0.09094145309269727 } },
      { { 0.000475993305781, 0.06923315329746019, 0.1345446587879024 } },
      { { 0.000526815118641, 0.09380751997100091, 0.1597986971590518 } },
      { { 0.000564304856051, 0.11959776925851553, 0.17623731960187097 } },
      { { 0.000591450107661, 0.1462944266839557, 0.18779859315577804 } },
      { { 0.000610456125787, 0.17372158857969944, 0.19638310871162132 } },
      { { 0.000623025286071, 0.20178171635925438, 0.20301334075208635 } },
      { { 0.000630561876176, 0.23042860761585268, 0.20828355715615282 } },
      { { 0.00063430927676, 0.25965307523608, 0.21255830404253653 } },
      { { 0.000517626894574, 0.08794707624386362, 0.05525153980013495 } },
      { { 0.000556484031331, 0.11204509589546444, 0.09159274175854999 } },
      { { 0.000585642667104, 0.13744522362543432, 0.11712413186869162 } },
      { { 0.000606638692578, 0.1638235897217632, 0.13599893913777458 } },
      { { 0.000620882496223, 0.19099329347026625, 0.1505076183425183 } },
      { { 0.000629631429782, 0.21885026631240181, 0.16199557887487917 } },
      { { 0.000634042375679, 0.24734548378169824, 0.17129282153427458 } },
      { { 0.000582962767711, 0.1327847126725995, 0.03988414867701234 } },
      { { 0.000604869337608, 0.1575226934924021, 0.06975403592073667 } },
      { { 0.000620236231773, 0.18331731309264804, 0.09286112723766339 } },
      { { 0.00062990053284, 0.20997811889533646, 0.11122182993593546 } },
      { { 0.000634772239061, 0.2373961100597402, 0.12612209724859638 } },
      { { 0.000620377898124, 0.17930414564418112, 0.03142567934766161 } },
      { { 0.000630841467124, 0.20435855567306369, 0.05670247592940609 } },
      { { 0.000636270646696, 0.23033583124044069, 0.07740125753569922 } },
      { { 0.000637541417033, 0.22665311278358263, 0.02611611008118614 } } }
  };

public:
  using Weights = Eigen::VectorXd;
  using Points = Eigen::Matrix3Xd;

  static auto quadrature() -> std::tuple<Weights, Points>
  {
    Weights ws = Weights::Zero(1730);
    Points xs = Points::Zero(3, 1730);

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
