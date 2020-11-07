#pragma once

#include "SpinOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

/** @class H_B_spin
 *
 * @brief Magnetic spin operator
 *
 * Interaction operator obtained by differentiating the spin Hamiltonian wrt
 * the external magnetic field B:
 *
 * dH/dB = H_B_dip + H_B_spin
 *
 * H_B_spin = -\sum_j m_j
 *
 * where m_j is the magnetic moment of the electron.
 */

class H_B_spin final : public RankOneTensorOperator<3> {
public:
    H_B_spin() {
        const double g_e = PHYSCONST::g_e;

        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &h = (*this);
        h[0] = (-g_e / 2.0) * s[0];
        h[1] = (-g_e / 2.0) * s[1];
        h[2] = (-g_e / 2.0) * s[2];
        h[0].name() = "h_B_spin[x]";
        h[1].name() = "h_B_spin[y]";
        h[2].name() = "h_B_spin[z]";
    }

private:
    SpinOperator s;
};

} // namespace mrchem
