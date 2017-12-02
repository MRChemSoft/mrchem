#pragma once

#include "QMTensorOperator.h"
#include "AngularMomentumOperator.h"

class H_B_dip : public RankOneTensorOperator<3> {
public:
    H_B_dip(DerivativeOperator<3> &d, const double *o = 0) : l(d, o) {
        initializeTensorOperator();
    }
    virtual ~H_B_dip() { }

protected:
    AngularMomentumOperator l;

    void initializeTensorOperator() {
        RankZeroTensorOperator &l_x = this->l[0];
        RankZeroTensorOperator &l_y = this->l[1];
        RankZeroTensorOperator &l_z = this->l[2];

        RankOneTensorOperator<3> &h = (*this);
        h[0] = 0.5*l_x;
        h[1] = 0.5*l_y;
        h[2] = 0.5*l_z;
    }
};


