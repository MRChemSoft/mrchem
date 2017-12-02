#pragma once

#include "QMTensorOperator.h"
#include "PositionOperator.h"
#include "MomentumOperator.h"

class AngularMomentumOperator : public RankOneTensorOperator<3> {
public:
    AngularMomentumOperator(DerivativeOperator<3> &D, const double *o = 0)
            : r(o), p(D) {
        initializeTensorOperator();
    }
    virtual ~AngularMomentumOperator() { }

protected:
    PositionOperator r;
    MomentumOperator p;

    void initializeTensorOperator() {
        RankZeroTensorOperator &r_x = this->r[0];
        RankZeroTensorOperator &r_y = this->r[1];
        RankZeroTensorOperator &r_z = this->r[2];
        RankZeroTensorOperator &p_x = this->p[0];
        RankZeroTensorOperator &p_y = this->p[1];
        RankZeroTensorOperator &p_z = this->p[2];

        RankOneTensorOperator<3> &h = (*this);
        h[0] = (r_y*p_z - r_z*p_y);
        h[1] = (r_z*p_x - r_x*p_z);
        h[2] = (r_x*p_y - r_y*p_x);
    }
};

