
#pragma once

#include <array>
#include <tuple>

#include <Eigen/Core>

#include "lebedev_utils.hpp"
#include "meta_utils.hpp"


namespace detail {
class Lebedev_077_2030 final
{
private:
  /** Number of points in A1 block */
  static constexpr auto N_a1 = 6;
  /** Quadrature weight in A1 block */
  static constexpr auto a1 = 4.6560318992e-05;

  /** Number of points in A3 block */
  static constexpr auto N_a3 = 8;
  /** Quadrature weight in A3 block */
  static constexpr auto a3 = 0.00054215491953;

  /** Number of points in each PQ0 block */
  static constexpr auto N_pq0 = 24;
  /** Quadrature weight(s) and point(s) in PQ0 block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as \f$(p, q,
   * 0) = (\sin(x * \pi), \cos(x * \pi), 0)\f$
   */
  static constexpr std::array<std::array<double, 2>, 6> pq0 = { { { { 0.000256800249773, 0.02382873604048719 } },
                                                                  { { 0.000382721170029, 0.05919454144498195 } },
                                                                  { { 0.000457949156192, 0.09897501240267485 } },
                                                                  { { 0.000504200396908, 0.14097739113409083 } },
                                                                  { { 0.000531270888998, 0.18418303972165878 } },
                                                                  { { 0.000543840179075, 0.22799849848688084 } } } };

  /** Number of points in each LLM block */
  static constexpr auto N_llm = 24;
  /** Quadrature weight(s) and point(s) in LLM block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(l, l, m) = (\frac{\sin(x * \pi)}{\sqrt{2}}, \frac{\sin(x * \pi)}{\sqrt{2}}, \cos(x * \pi))\f$
   */
  static constexpr std::array<std::array<double, 2>, 18> llm = { { { { 0.000177852213335, 0.011440240333553339 } },
                                                                   { { 0.000281132540568, 0.02884654222301814 } },
                                                                   { { 0.000354889631263, 0.04918482752065367 } },
                                                                   { { 0.000409031089717, 0.07129976191437562 } },
                                                                   { { 0.000449328613417, 0.09462257239966308 } },
                                                                   { { 0.000479372844796, 0.11883307451319587 } },
                                                                   { { 0.000501541531916, 0.14374040519867876 } },
                                                                   { { 0.000517512737268, 0.16922931750291914 } },
                                                                   { { 0.000528552226208, 0.19523304655136767 } },
                                                                   { { 0.000535683270371, 0.2217185876480012 } },
                                                                   { { 0.000539791473618, 0.24867829648769377 } },
                                                                   { { 0.00054168994416, 0.2761247273453157 } },
                                                                   { { 0.000541930847689, 0.33260491862951663 } },
                                                                   { { 0.000541693690203, 0.361724274183216 } },
                                                                   { { 0.00054195443387, 0.39148038327329115 } },
                                                                   { { 0.000542898365663, 0.42187592637813365 } },
                                                                   { { 0.00054422865001, 0.452847070422067 } },
                                                                   { { 0.000545225034506, 0.48423160979755336 } } } };

  /** Number of points in each RSW block */
  static constexpr auto N_rsw = 48;
  /** Quadrature weight(s) and point(s) in RSW block.
   *
   * Each element is a {weight, x, y} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(r, s, w) = (\sin(x * \pi)\cos(y * \pi), \sin(x * \pi)\sin(y * \pi), \cos(x*\pi))\f$
   */
  static constexpr std::array<std::array<double, 3>, 30> rsw = {
    { { { 0.00033160418732, 0.041717375566021495, 0.0909121294226695 } },
      { { 0.000389911356715, 0.0623698186564717, 0.1344754439526608 } },
      { { 0.00043433433272, 0.08470338444769895, 0.1596897782776591 } },
      { { 0.000467941526232, 0.1081870685924007, 0.1760888376667163 } },
      { { 0.000493084798163, 0.1325248525263869, 0.1876123327182394 } },
      { { 0.000511503186754, 0.1575419782291847, 0.19616323676607925 } },
      { { 0.000524521714846, 0.18313463768637675, 0.2027669357714648 } },
      { { 0.00053320414999, 0.20924485380242722, 0.20802142331344914 } },
      { { 0.000538458312602, 0.23584693931713388, 0.21229618706426556 } },
      { { 0.00054110672108, 0.2629397299261896, 0.21582896181460065 } },
      { { 0.000425979739147, 0.07939622760508908, 0.055182374747974626 } },
      { { 0.000460493136846, 0.10134657952582082, 0.09144055632471736 } },
      { { 0.000487181487826, 0.1245144128130184, 0.1168857298537963 } },
      { { 0.000507224291007, 0.1485925561541061, 0.13567815467589486 } },
      { { 0.000521706984524, 0.1733972541505497, 0.15011458258626764 } },
      { { 0.000531578596628, 0.1988190002039459, 0.16154737711489003 } },
      { { 0.000537683370876, 0.22479687342928137, 0.17081525014445756 } },
      { { 0.000540803209207, 0.2513043794619103, 0.17845916790513258 } },
      { { 0.00048427449179, 0.12031244726303181, 0.03977490130136878 } },
      { { 0.000504892607619, 0.14292449218331743, 0.06951809343561059 } },
      { { 0.000520260798048, 0.16651144408744856, 0.09249764235304189 } },
      { { 0.000530993238833, 0.19088638537196473, 0.11074376794678162 } },
      { { 0.00053774197709, 0.2159349509430741, 0.12555631421581928 } },
      { { 0.000541169633168, 0.2415923854202134, 0.1377912289972581 } },
      { { 0.000519799629328, 0.1629280674120531, 0.031277144226667615 } },
      { { 0.000531112083662, 0.18589222528429114, 0.056391299211738244 } },
      { { 0.000538430931996, 0.20969245446347537, 0.07693651699435375 } },
      { { 0.000542185950405, 0.23421699784124456, 0.09400845965175704 } },
      { { 0.000539094835505, 0.20645400832343744, 0.02593526986849543 } },
      { { 0.000543331270503, 0.22958259580795992, 0.04771629680091758 } } }
  };

public:
  using Weights = Eigen::VectorXd;
  using Points = Eigen::Matrix3Xd;

  static auto quadrature() -> std::tuple<Weights, Points>
  {
    Weights ws = Weights::Zero(2030);
    Points xs = Points::Zero(3, 2030);

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

