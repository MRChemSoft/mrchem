
#pragma once

#include <array>
#include <tuple>

#include <Eigen/Core>

#include "lebedev_utils.hpp"
#include "meta_utils.hpp"


namespace detail {
class Lebedev_083_2354 final
{
private:
  /** Number of points in A1 block */
  static constexpr auto N_a1 = 6;
  /** Quadrature weight in A1 block */
  static constexpr auto a1 = 3.9226162707e-05;

  /** Number of points in A2 block */
  static constexpr auto N_a2 = 12;
  /** Quadrature weight in A2 block */
  static constexpr auto a2 = 0.000470383175085;

  /** Number of points in A3 block */
  static constexpr auto N_a3 = 8;
  /** Quadrature weight in A3 block */
  static constexpr auto a3 = 0.000467820280128;

  /** Number of points in each PQ0 block */
  static constexpr auto N_pq0 = 24;
  /** Quadrature weight(s) and point(s) in PQ0 block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as \f$(p, q,
   * 0) = (\sin(x * \pi), \cos(x * \pi), 0)\f$
   */
  static constexpr std::array<std::array<double, 2>, 6> pq0 = { { { { 0.000209994228107, 0.02148323008071945 } },
                                                                  { { 0.000317226915071, 0.053688384512778804 } },
                                                                  { { 0.000383205135855, 0.09012092708102555 } },
                                                                  { { 0.000425219381815, 0.12874501337985692 } },
                                                                  { { 0.000451380796375, 0.16860448658816918 } },
                                                                  { { 0.000465779746911, 0.209148252869559 } } } };

  /** Number of points in each LLM block */
  static constexpr auto N_llm = 24;
  /** Quadrature weight(s) and point(s) in LLM block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(l, l, m) = (\frac{\sin(x * \pi)}{\sqrt{2}}, \frac{\sin(x * \pi)}{\sqrt{2}}, \cos(x * \pi))\f$
   */
  static constexpr std::array<std::array<double, 2>, 19> llm = { { { { 0.000143783222898, 0.0103105356590626 } },
                                                                   { { 0.000230357249358, 0.02604407895727144 } },
                                                                   { { 0.000293311075245, 0.04454481490925946 } },
                                                                   { { 0.000340290599836, 0.06473491178346792 } },
                                                                   { { 0.000375913846687, 0.08607922577105943 } },
                                                                   { { 0.00040306384479, 0.10827224463759817 } },
                                                                   { { 0.000423659143224, 0.13112750294976167 } },
                                                                   { { 0.000439052265695, 0.15452812760064788 } },
                                                                   { { 0.000450252346663, 0.17840164434395828 } },
                                                                   { { 0.000458057772778, 0.20270607993191486 } },
                                                                   { { 0.000463139161662, 0.22742189873716076 } },
                                                                   { { 0.00046609289537, 0.2525471144673109 } },
                                                                   { { 0.000467475180794, 0.278094061578268 } },
                                                                   { { 0.000467641490393, 0.33055750758729435 } },
                                                                   { { 0.000467408649235, 0.35754201306062255 } },
                                                                   { { 0.000467492853948, 0.38506984706068115 } },
                                                                   { { 0.000468074897969, 0.413149299426279 } },
                                                                   { { 0.000469044980639, 0.4417452451509824 } },
                                                                   { { 0.000469987707586, 0.4707553754870439 } } } };

  /** Number of points in each RSW block */
  static constexpr auto N_rsw = 48;
  /** Quadrature weight(s) and point(s) in RSW block.
   *
   * Each element is a {weight, x, y} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(r, s, w) = (\sin(x * \pi)\cos(y * \pi), \sin(x * \pi)\sin(y * \pi), \cos(x*\pi))\f$
   */
  static constexpr std::array<std::array<double, 3>, 36> rsw = {
    { { { 0.000273336280052, 0.037743625637352785, 0.09088472512920288 } },
      { { 0.000323548536846, 0.05658124612192914, 0.13442239982419044 } },
      { { 0.000362490872601, 0.07700624085599411, 0.15960668222432448 } },
      { { 0.000392554007071, 0.09852377819321616, 0.17597464333125812 } },
      { { 0.000415612978112, 0.12085224031518965, 0.18746763179651368 } },
      { { 0.000433064498462, 0.14382072327200804, 0.19599015346591778 } },
      { { 0.000445967772592, 0.16732318186207223, 0.20256938328271545 } },
      { { 0.000455159300446, 0.19129519268422873, 0.2078055112217531 } },
      { { 0.000461334146275, 0.21570120393595657, 0.2120708126146701 } },
      { { 0.000465101961827, 0.24052715022858334, 0.2156066524708237 } },
      { { 0.00046702495361, 0.26577591697740255, 0.2185744280602722 } },
      { { 0.000354955557644, 0.07216504764039339, 0.055130070446995944 } },
      { { 0.000385610824525, 0.09228070552880042, 0.0913242194519082 } },
      { { 0.000409862284576, 0.11354199196436619, 0.11670182943449679 } },
      { { 0.000428632860427, 0.1356585148603882, 0.13542796525954431 } },
      { { 0.000442780219899, 0.15845204047317685, 0.14980340556108718 } },
      { { 0.000453047351149, 0.18181171606183896, 0.16118476007373866 } },
      { { 0.00046008054757, 0.20567040966642677, 0.1704157879330728 } },
      { { 0.000464459905996, 0.2299914766537107, 0.17804384755163255 } },
      { { 0.000466727445571, 0.2547609952738979, 0.1844353087637181 } },
      { { 0.000406936051802, 0.10971921112740224, 0.03969101287175254 } },
      { { 0.000426044281992, 0.1305095211763083, 0.06933542701583828 } },
      { { 0.000440867850803, 0.15220870821515906, 0.0922129333437181 } },
      { { 0.000451874811555, 0.17463595623685488, 0.1103629057453786 } },
      { { 0.000459556487538, 0.19767599129304933, 0.12509380222310257 } },
      { { 0.000464398877432, 0.22125786065218556, 0.13727104621488875 } },
      { { 0.000466882749165, 0.2453423661150161, 0.14748176427821896 } },
      { { 0.000440054182374, 0.1489712500225431, 0.031161557211991396 } },
      { { 0.000451451289019, 0.1701392349394957, 0.05614611735613247 } },
      { { 0.000459619862735, 0.1920763306362031, 0.07656320836273045 } },
      { { 0.00046486590168, 0.21466903585235447, 0.0935232678749094 } },
      { { 0.000467550201716, 0.237846060277366, 0.10780208571690666 } },
      { { 0.000459849447646, 0.18918326318817155, 0.025791431603201697 } },
      { { 0.000465491695515, 0.21055131325424792, 0.04742096375190213 } },
      { { 0.000468470977951, 0.23260332842524298, 0.06577176881442136 } },
      { { 0.000469144553911, 0.2298980098880378, 0.022118446350596878 } } }
  };

public:
  using Weights = Eigen::VectorXd;
  using Points = Eigen::Matrix3Xd;

  static auto quadrature() -> std::tuple<Weights, Points>
  {
    Weights ws = Weights::Zero(2354);
    Points xs = Points::Zero(3, 2354);

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

