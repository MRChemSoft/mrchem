#pragma once

#include "qmoperators/RankTwoTensorOperator.h"
#include "PositionOperator.h"

namespace mrchem {

class H_BB_dia final : public RankTwoTensorOperator<3,3> {
public:
    H_BB_dia(const double *o = 0)
            : r(o) {
        RankTwoTensorOperator<3,3> &h = (*this);
        h[0][0] = ( 1.0/4.0)*(r[1]*r[1] + r[2]*r[2]);
        h[0][1] = (-1.0/4.0)*(r[0]*r[1]);
        h[0][2] = (-1.0/4.0)*(r[0]*r[2]);
        h[1][0] = (-1.0/4.0)*(r[1]*r[0]);
        h[1][1] = ( 1.0/4.0)*(r[0]*r[0] + r[2]*r[2]);
        h[1][2] = (-1.0/4.0)*(r[1]*r[2]);
        h[2][0] = (-1.0/4.0)*(r[2]*r[0]);
        h[2][1] = (-1.0/4.0)*(r[2]*r[1]);
        h[2][2] = ( 1.0/4.0)*(r[0]*r[0] + r[1]*r[1]);
    }
    ~H_BB_dia() { }

protected:
    PositionOperator r;
};

} //namespace mrchem
