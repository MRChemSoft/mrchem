#pragma once

#include "RankOneTensorOperator.h"
#include "PositionOperator.h"
#include "DistanceOperator.h"
#include "NuclearOperator.h"

namespace mrchem {

class X_rm3 final : public RankOneTensorOperator<3> {
public:
    X_rm3(const Nuclei &nucs, const double *R_k = 0)
	: r_m3(3.0, R_k, 1.0e-3),
	    r(nucs[0].getCoord()) {
        RankOneTensorOperator<3> &h = (*this);
        h[0] = r_m3*r[0];
        h[1] = r_m3*r[1];
        h[2] = r_m3*r[2];
    }
    ~X_rm3() { }

protected:
    DistanceOperator r_m3;
    PositionOperator r;
};

} //namespace mrchem
