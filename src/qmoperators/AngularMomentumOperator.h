#pragma once

#include "PositionOperator.h"
#include "MomentumOperator.h"

namespace mrchem {

class AngularMomentumOperator final : public RankOneTensorOperator<3> {
public:
    AngularMomentumOperator(mrcpp::DerivativeOperator<3> &D, const double *o = 0)
            : r(o),
              p(D) {
        RankOneTensorOperator<3> &h = (*this);
        h[0] = (r[1]*p[2] - r[2]*p[1]);
        h[1] = (r[2]*p[0] - r[0]*p[2]);
        h[2] = (r[0]*p[1] - r[1]*p[0]);
    }
    ~AngularMomentumOperator() { }

protected:
    PositionOperator r;
    MomentumOperator p;
};

} //namespace mrchem
