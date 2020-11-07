#pragma once

#include "AngularMomentumOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

/** @class H_B_dip
 *
 * @brief Magnetic dipole operator
 *
 * Interaction operator obtained by differentiating the spin Hamiltonian wrt
 * the external magnetic field B:
 *
 * dH/dB = H_B_dip + H_B_spin
 *
 * H_B_dip = \frac{1}{2} \sum_j l_{jO}
 *
 * where l_{jO} is the orbital angular momentum.
 */

class H_B_dip final : public RankOneTensorOperator<3> {
public:
    H_B_dip(std::shared_ptr<mrcpp::DerivativeOperator<3>> D, const mrcpp::Coord<3> &o)
            : l(D, o) {
        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &h = (*this);
        h[0] = 0.5 * l[0];
        h[1] = 0.5 * l[1];
        h[2] = 0.5 * l[2];
        h[0].name() = "h_B_dip[x]";
        h[1].name() = "h_B_dip[y]";
        h[2].name() = "h_B_dip[z]";
    }

private:
    AngularMomentumOperator l;
};

} // namespace mrchem
