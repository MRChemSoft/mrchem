#pragma once

#include "RankOneTensorOperator.h"
#include "SpinOperator.h"

namespace mrchem {

class H_B_spin final : public RankOneTensorOperator<3> {
public:
    H_B_spin() {
        const double g_e = PHYSCONST::g_e;

        RankOneTensorOperator<3> &h = (*this);
        h[0] = (-g_e/2.0)*s[0];
        h[1] = (-g_e/2.0)*s[1];
        h[2] = (-g_e/2.0)*s[2];
    }
    ~H_B_spin() { }

protected:
    SpinOperator s;
};

} //namespace mrchem
