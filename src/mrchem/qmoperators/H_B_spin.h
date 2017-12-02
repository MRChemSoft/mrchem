#pragma once

#include "QMTensorOperator.h"
#include "SpinOperator.h"

class H_B_spin : public RankOneTensorOperator<3> {
public:
    H_B_spin() { initializeTensorOperator(); }
    virtual ~H_B_spin() { }
protected:
    SpinOperator s;

    void initializeTensorOperator() {
        static double g_e = 2.00231930436182;
        RankZeroTensorOperator &s_x = this->s[0];
        RankZeroTensorOperator &s_y = this->s[1];
        RankZeroTensorOperator &s_z = this->s[2];

        RankOneTensorOperator<3> &h = (*this);
        h[0] = (-g_e/2.0)*s_x;
        h[1] = (-g_e/2.0)*s_y;
        h[2] = (-g_e/2.0)*s_z;
    }
};


