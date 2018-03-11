#pragma once

#include "RankOneTensorOperator.h"
#include "AngularMomentumOperator.h"

namespace mrchem {

class H_B_dip final : public RankOneTensorOperator<3> {
public:
    H_B_dip(mrcpp::DerivativeOperator<3> &D, const double *o = 0)
            : l(D, o) {
        RankOneTensorOperator<3> &h = (*this);
        h[0] = 0.5*l[0];
        h[1] = 0.5*l[1];
        h[2] = 0.5*l[2];
    }
    ~H_B_dip() { }

protected:
    AngularMomentumOperator l;
};

} //namespace mrchem
