#pragma once

#include "QMTensorOperator.h"
#include "PositionOperator.h"

class H_BB_dia : public RankTwoTensorOperator<3,3> {
public:
    H_BB_dia(const double *o = 0) : r(o) {
        initializeTensorOperator();
    }
    virtual ~H_BB_dia() { }

protected:
    PositionOperator r;

    void initializeTensorOperator() {
        RankZeroTensorOperator &r_x = this->r[0];
        RankZeroTensorOperator &r_y = this->r[1];
        RankZeroTensorOperator &r_z = this->r[2];

        RankTwoTensorOperator<3,3> &h = (*this);
        h[0][0] =  1.0/4.0*(r_y*r_y + r_z*r_z);
        h[0][1] = -1.0/4.0*(r_x*r_y);
        h[0][2] = -1.0/4.0*(r_x*r_z);
        h[1][0] = -1.0/4.0*(r_y*r_x);
        h[1][1] =  1.0/4.0*(r_x*r_x + r_z*r_z);
        h[1][2] = -1.0/4.0*(r_y*r_z);
        h[2][0] = -1.0/4.0*(r_z*r_x);
        h[2][1] = -1.0/4.0*(r_z*r_y);
        h[2][2] =  1.0/4.0*(r_x*r_x + r_y*r_y);
    }
};


