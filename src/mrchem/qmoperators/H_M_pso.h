#pragma once

#include "QMTensorOperator.h"
#include "AngularMomentumOperator.h"
#include "CubicPotential.h"

class H_M_pso : public RankOneTensorOperator<3> {
public:
    H_M_pso(DerivativeOperator<3> &d, const double *k = 0) : r_m3(1.0, k), l(d, k) {
        initializeTensorOperator();
    }
    virtual ~H_M_pso() { }

protected:
    CubicPotential r_m3;
    AngularMomentumOperator l;

    void initializeTensorOperator() {
        static double alpha = 7.2973525664;
        RankZeroTensorOperator &l_x = this->l[0];
        RankZeroTensorOperator &l_y = this->l[1];
        RankZeroTensorOperator &l_z = this->l[2];

        RankOneTensorOperator<3> &h = (*this);
        h[0] = (alpha*alpha)*r_m3*l_x;
        h[1] = (alpha*alpha)*r_m3*l_y;
        h[2] = (alpha*alpha)*r_m3*l_z;
    }
};


