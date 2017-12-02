#pragma once

#include "QMTensorOperator.h"
#include "DeltaOperator.h"
#include "H_B_spin.h"

class H_M_fc : public RankOneTensorOperator<3> {
public:
    H_M_fc(const double *o = 0) : delta(o, 1.0e6) {
        initializeTensorOperator();
    }
    virtual ~H_M_fc() { }

protected:
    H_B_spin s;
    DeltaOperator delta;

    void initializeTensorOperator() {
        static double coef = -8.0*pi/3.0;
        static double alpha = 7.2973525664;
        RankZeroTensorOperator &s_x = this->s[0];
        RankZeroTensorOperator &s_y = this->s[1];
        RankZeroTensorOperator &s_z = this->s[2];

        RankOneTensorOperator<3> &h = (*this);
        h[0] = (coef*alpha*alpha)*delta*s_x;
        h[1] = (coef*alpha*alpha)*delta*s_y;
        h[2] = (coef*alpha*alpha)*delta*s_z;
    }
};


