#pragma once

#include "RankOneTensorOperator.h"
#include "AngularMomentumOperator.h"
#include "DistanceOperator.h"

namespace mrchem {

class H_M_pso final : public RankOneTensorOperator<3> {
public:
    H_M_pso(mrcpp::DerivativeOperator<3> &D, const double *k = 0)
            : r_m3(3.0, k),
              l(D, k) {
        const double alpha_2 = PHYSCONST::alpha * PHYSCONST::alpha;
        RankOneTensorOperator<3> &h = (*this);
        h[0] = (alpha_2)*r_m3*l[0];
        h[1] = (alpha_2)*r_m3*l[1];
        h[2] = (alpha_2)*r_m3*l[2];
    }
    ~H_M_pso() { }

protected:
    DistanceOperator r_m3;
    AngularMomentumOperator l;
};

} //namespace mrchem
