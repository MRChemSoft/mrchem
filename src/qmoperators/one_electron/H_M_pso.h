#pragma once

#include "MomentumOperator.h"
#include "NuclearGradientOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

/** @class H_M_pso
 *
 * @brief Paramagnetic Spin-Orbit operator
 *
 * Interaction operator obtained by differentiating the spin Hamiltonian wrt
 * the nuclear magnetic moment of nucleaus K:
 *
 * dH/dM_K = H_M_pso + H_M_sd + H_M_fc
 *
 * H_M_pso = alpha^2 \sum_j \frac{l_{jK}}{r_{jK}^3}
 *
 * where l_{jK} is the orbital angular momentum.
 */

class H_M_pso final : public RankOneTensorOperator<3> {
public:
    H_M_pso(std::shared_ptr<mrcpp::DerivativeOperator<3>> D, const mrcpp::Coord<3> &k, double c)
            : p(D)
            , r_m3(1.0, k, c) {
        const double alpha_2 = PHYSCONST::alpha * PHYSCONST::alpha;

        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &h = (*this);
        h[0] = alpha_2 * (r_m3[1] * p[2] - r_m3[2] * p[1]);
        h[1] = alpha_2 * (r_m3[2] * p[0] - r_m3[0] * p[2]);
        h[2] = alpha_2 * (r_m3[0] * p[1] - r_m3[1] * p[0]);
        h[0].name() = "h_M_pso[x]";
        h[1].name() = "h_M_pso[y]";
        h[2].name() = "h_M_pso[z]";
    }

private:
    MomentumOperator p;
    NuclearGradientOperator r_m3;
};

} // namespace mrchem
