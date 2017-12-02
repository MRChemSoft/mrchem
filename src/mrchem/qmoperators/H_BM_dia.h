#pragma once

#include "QMTensorOperator.h"
#include "PositionOperator.h"
#include "CubicPotential.h"

class H_BM_dia : public RankTwoTensorOperator<3,3> {
public:
    H_BM_dia(const double *o = 0, const double *k = 0)
            : r_o(o), r_k(k), r_m3(1.0, k) {
        initializeTensorOperator();
    }
    virtual ~H_BM_dia() { }

protected:
    PositionOperator r_o;
    PositionOperator r_k;
    CubicPotential r_m3;

    void initializeTensorOperator() {
        static double alpha = 7.2973525664;

        RankZeroTensorOperator &o_x = this->r_o[0];
        RankZeroTensorOperator &o_y = this->r_o[1];
        RankZeroTensorOperator &o_z = this->r_o[2];
        RankZeroTensorOperator &k_x = this->r_k[0];
        RankZeroTensorOperator &k_y = this->r_k[1];
        RankZeroTensorOperator &k_z = this->r_k[2];

        RankTwoTensorOperator<3,3> &h = (*this);
        h[0][0] = -(alpha*alpha/2.0)*r_m3*(o_y*k_y + o_z*k_z);
        h[0][1] =  (alpha*alpha/2.0)*r_m3*(o_x*k_y);
        h[0][2] =  (alpha*alpha/2.0)*r_m3*(o_x*k_z);
        h[1][0] =  (alpha*alpha/2.0)*r_m3*(o_y*k_x);
        h[1][1] = -(alpha*alpha/2.0)*r_m3*(o_x*k_x + o_z*k_z);
        h[1][2] =  (alpha*alpha/2.0)*r_m3*(o_y*k_z);
        h[2][0] =  (alpha*alpha/2.0)*r_m3*(o_z*k_x);
        h[2][1] =  (alpha*alpha/2.0)*r_m3*(o_z*k_y);
        h[2][2] = -(alpha*alpha/2.0)*r_m3*(o_x*k_x + o_y*k_y);
    }
};


