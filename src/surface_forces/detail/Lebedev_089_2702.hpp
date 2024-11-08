
#pragma once

#include <array>
#include <tuple>

#include <Eigen/Core>

#include "lebedev_utils.hpp"
#include "meta_utils.hpp"


namespace detail {
class Lebedev_089_2702 final
{
private:
  /** Number of points in A1 block */
  static constexpr auto N_a1 = 6;
  /** Quadrature weight in A1 block */
  static constexpr auto a1 = 2.9986751499e-05;

  /** Number of points in A3 block */
  static constexpr auto N_a3 = 8;
  /** Quadrature weight in A3 block */
  static constexpr auto a3 = 0.00040778605295;

  /** Number of points in each PQ0 block */
  static constexpr auto N_pq0 = 24;
  /** Quadrature weight(s) and point(s) in PQ0 block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as \f$(p, q,
   * 0) = (\sin(x * \pi), \cos(x * \pi), 0)\f$
   */
  static constexpr std::array<std::array<double, 2>, 7> pq0 = { { { { 0.000173898681175, 0.019504083272014666 } },
                                                                  { { 0.000265961604528, 0.048992350041236 } },
                                                                  { { 0.000324059600817, 0.0825332235620925 } },
                                                                  { { 0.000362119596443, 0.11822519839431116 } },
                                                                  { { 0.000386883833076, 0.15516757030789682 } },
                                                                  { { 0.000401891153269, 0.19284298215381002 } },
                                                                  { { 0.000408992943298, 0.23090699168193182 } } } };

  /** Number of points in each LLM block */
  static constexpr auto N_llm = 24;
  /** Quadrature weight(s) and point(s) in LLM block.
   *
   * Each element is a {weight, x} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(l, l, m) = (\frac{\sin(x * \pi)}{\sqrt{2}}, \frac{\sin(x * \pi)}{\sqrt{2}}, \cos(x * \pi))\f$
   */
  static constexpr std::array<std::array<double, 2>, 21> llm = {
    { { { 0.000118534919252, 0.009299621174798444 } }, { { 0.000191340864343, 0.023659215108291817 } },
      { { 0.000245288657721, 0.040594837447854165 } }, { { 0.000286240818329, 0.05913124047454959 } },
      { { 0.000317803225826, 0.07877165755876729 } },  { { 0.000342294566763, 0.09922627578196008 } },
      { { 0.000361279052024, 0.1203145040693523 } },   { { 0.000375863822982, 0.14192006474436658 } },
      { { 0.000386871179886, 0.16396770467049449 } },  { { 0.000394942993319, 0.18641014672799455 } },
      { { 0.000400606810754, 0.20922039460299993 } },  { { 0.000404319214967, 0.23238702174833858 } },
      { { 0.000406494749581, 0.255911158481507 } },    { { 0.000407524561981, 0.27980436335312925 } },
      { { 0.000407642354089, 0.3287844859640549 } },   { { 0.000407428086225, 0.3539262986287412 } },
      { { 0.000407416375601, 0.37953686493412986 } },  { { 0.000407764779507, 0.4056266456277628 } },
      { { 0.000408451755278, 0.4321769981995628 } },   { { 0.000409246845922, 0.4591229832325327 } },
      { { 0.000409787268724, 0.4863416665034887 } } }
  };

  /** Number of points in each RSW block */
  static constexpr auto N_rsw = 48;
  /** Quadrature weight(s) and point(s) in RSW block.
   *
   * Each element is a {weight, x, y} tuple. The actual quadrature point in Cartesian coordinates is computed as
   * \f$(r, s, w) = (\sin(x * \pi)\cos(y * \pi), \sin(x * \pi)\sin(y * \pi), \cos(x*\pi))\f$
   */
  static constexpr std::array<std::array<double, 3>, 42> rsw = {
    { { { 0.000227990752771, 0.03436681470127729, 0.09086849544871468 } },
      { { 0.000271520549058, 0.051644421330199426, 0.13438217012923867 } },
      { { 0.00030579178967, 0.07042550801173678, 0.15954236093630902 } },
      { { 0.000332691305245, 0.09024803669564561, 0.17588550976799394 } },
      { { 0.000353733471189, 0.11084459756408502, 0.18735370470068727 } },
      { { 0.000370056750078, 0.13204947781981466, 0.19585248971516858 } },
      { { 0.000382524537259, 0.1537567152755291, 0.2024101946225456 } },
      { { 0.000391812517152, 0.1758985905010745, 0.20762836888005923 } },
      { { 0.000398472041994, 0.1984336588948514, 0.21188095922911135 } },
      { { 0.000402974600334, 0.22133971967495242, 0.21541144144517566 } },
      { { 0.000405742863216, 0.24460950296105605, 0.21838392518490038 } },
      { { 0.000407171927411, 0.2682478541828218, 0.22091169113093448 } },
      { { 0.000299023695066, 0.06598188689594922, 0.055089534111938475 } },
      { { 0.000326295173421, 0.08451357920461729, 0.09123376892849176 } },
      { { 0.000348263460824, 0.10412844646981298, 0.11655782680613067 } },
      { { 0.00036565966817, 0.12455231517854437, 0.13523033141502266 } },
      { { 0.000379174046779, 0.14561357152343674, 0.14955486971484258 } },
      { { 0.000389403445016, 0.16720229214098942, 0.16089079428781922 } },
      { { 0.000396860024551, 0.18924841946053073, 0.17008496585779787 } },
      { { 0.000401993135142, 0.2117094049587665, 0.17768852416441078 } },
      { { 0.000405210880128, 0.2345628880915592, 0.1840726291968587 } },
      { { 0.000406897861394, 0.2578021611634835, 0.1894935735877868 } },
      { { 0.000345427535132, 0.10062502552805537, 0.0396255830809811 } },
      { { 0.000362996353701, 0.11983775907421261, 0.06919193291928188 } },
      { { 0.000377018723389, 0.13990431349250007, 0.09198727485723154 } },
      { { 0.000387860861369, 0.16065107710970586, 0.11005739257320142 } },
      { { 0.000395906527022, 0.18196424996532481, 0.12471640471346784 } },
      { { 0.000401528697546, 0.20377022223408525, 0.1368356052613229 } },
      { { 0.000405086678561, 0.22602383303298287, 0.1470089338642573 } },
      { { 0.000406932018505, 0.24870114433171436, 0.1556502756669482 } },
      { { 0.000376012096406, 0.13695207931460154, 0.03107036643743837 } },
      { { 0.000387096956442, 0.15656068154418076, 0.05595093547248068 } },
      { { 0.000395528779053, 0.17688469742993385, 0.07626210978060294 } },
      { { 0.00040153619113, 0.19781170620276076, 0.09312413288736685 } },
      { { 0.000405383698672, 0.2192674038421525, 0.10732270843070195 } },
      { { 0.00040735786733, 0.2412055444632103, 0.11941677099067832 } },
      { { 0.000395462837923, 0.17426866411547765, 0.025676175296736718 } },
      { { 0.000401764550885, 0.19410348608945335, 0.04718083967949124 } },
      { { 0.000405903034865, 0.21456648708853332, 0.06541131103330632 } },
      { { 0.000408056580948, 0.23558660319741667, 0.08102731647185069 } },
      { { 0.000406301875366, 0.2121495916641633, 0.02198403859038884 } },
      { { 0.00040871912928, 0.23209355666287107, 0.04097218051583061 } } }
  };

public:
  using Weights = Eigen::VectorXd;
  using Points = Eigen::Matrix3Xd;

  static auto quadrature() -> std::tuple<Weights, Points>
  {
    Weights ws = Weights::Zero(2702);
    Points xs = Points::Zero(3, 2702);

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

