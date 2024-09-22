
#pragma once

#include <array>
#include <tuple>

#include <Eigen/Core>

#include "lebedev_utils.hpp"
#include "meta_utils.hpp"

namespace detail {
class Lebedev_065_1454 final
{
private:
  /** Number of points in A1 block */
  static constexpr auto N_a1 = 6;
  /** Quadrature weight in A1 block */
  static constexpr auto a1 = 7.7771607433e-05;

  /** Number of points in A3 block */
  static constexpr auto N_a3 = 8;
  /** Quadrature weight in A3 block */
  static constexpr auto a3 = 0.0007557646413;

  /** Number of points in each PQ0 block */
  static constexpr auto N_pq0 = 24;
  /** Quadrature weight(s) and point(s) in PQ0 block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as \f$(p, q,
   * 0) = (\sin(x * \pi), \cos(x * \pi), 0)\f$
   */
  static constexpr std::array<std::array<double, 2>, 5> pq0 = { { { { 0.000402168044787, 0.03002483397375312 } },
                                                                  { { 0.000580487179395, 0.07356619587204107 } },
                                                                  { { 0.000679215195595, 0.1218942218184075 } },
                                                                  { { 0.000733674121129, 0.17244770702224507 } },
                                                                  { { 0.000758186630099, 0.22404572432210743 } } } };

  /** Number of points in each LLM block */
  static constexpr auto N_llm = 24;
  /** Quadrature weight(s) and point(s) in LLM block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(l, l, m) = (\frac{\sin(x * \pi)}{\sqrt{2}}, \frac{\sin(x * \pi)}{\sqrt{2}}, \cos(x * \pi))\f$
   */
  static constexpr std::array<std::array<double, 2>, 15> llm = { { { { 0.000284163380609, 0.014541973303948906 } },
                                                                   { { 0.000437441912705, 0.036256356902763215 } },
                                                                   { { 0.000541717474087, 0.06134347429283711 } },
                                                                   { { 0.000614800089136, 0.08841663566770758 } },
                                                                   { { 0.00066643944858, 0.11682578562312121 } },
                                                                   { { 0.000702503935692, 0.14622485492090825 } },
                                                                   { { 0.000726851178925, 0.176424040909554 } },
                                                                   { { 0.000742263753421, 0.20732548358335812 } },
                                                                   { { 0.000750954503584, 0.23889218046843466 } },
                                                                   { { 0.000754853505772, 0.27113178450101194 } },
                                                                   { { 0.000755408896977, 0.3378245651771407 } },
                                                                   { { 0.000755314717444, 0.3724214055050285 } },
                                                                   { { 0.000756476765329, 0.4079277989979188 } },
                                                                   { { 0.000758799180852, 0.4443064433703368 } },
                                                                   { { 0.000760826183203, 0.48135058138208653 } } } };

  /** Number of points in each RSW block */
  static constexpr auto N_rsw = 48;
  /** Quadrature weight(s) and point(s) in RSW block.
   *
   * Each element is a {weight, x, y} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(r, s, w) = (\sin(x * \pi)\cos(y * \pi), \sin(x * \pi)\sin(y * \pi), \cos(x*\pi))\f$
   */
  static constexpr std::array<std::array<double, 3>, 20> rsw = {
    { { { 0.00075382578598, 0.2558003577809087, 0.20859614189791004 } },
      { { 0.000748351724705, 0.22408472579890718, 0.2033204671768711 } },
      { { 0.000737176366111, 0.19306103334458685, 0.19666471662424795 } },
      { { 0.000718344889576, 0.16275432818948848, 0.18804149522722452 } },
      { { 0.000689581552982, 0.13324870864184388, 0.17643341386541955 } },
      { { 0.000648010580179, 0.10471633615407486, 0.15994385152980334 } },
      { { 0.000589755889659, 0.07747808165091459, 0.1346384497079675 } },
      { { 0.000509570884925, 0.05214524403629088, 0.09098925491262612 } },
      { { 0.000753690642891, 0.24277387514946042, 0.16254595611164124 } },
      { { 0.000747250596558, 0.21198920799674287, 0.15100683461207687 } },
      { { 0.000734301713228, 0.18199417202541565, 0.1364152050357355 } },
      { { 0.000713087158218, 0.1528797709057363, 0.1174379770106918 } },
      { { 0.000681702203211, 0.12483119354960753, 0.09179496313058289 } },
      { { 0.00063809411456, 0.09818737315696387, 0.055343676032884016 } },
      { { 0.000755038137792, 0.23269927807800983, 0.11182321068294632 } },
      { { 0.000747864664014, 0.203313866757253, 0.09333015692137783 } },
      { { 0.00073359187206, 0.17490196727498353, 0.07006345916885712 } },
      { { 0.000711012052766, 0.1476521039744757, 0.04002895490431486 } },
      { { 0.000757136397869, 0.22626624665276363, 0.057100620618486995 } },
      { { 0.000748990832908, 0.19874998248857403, 0.03161917277066387 } } }
  };

public:
  using Weights = Eigen::VectorXd;
  using Points = Eigen::Matrix3Xd;

  static auto quadrature() -> std::tuple<Weights, Points>
  {
    Weights ws = Weights::Zero(1454);
    Points xs = Points::Zero(3, 1454);

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
