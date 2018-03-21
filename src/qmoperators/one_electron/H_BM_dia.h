#pragma once

#include "RankTwoTensorOperator.h"
#include "PositionOperator.h"
#include "PowerPotential.h"

namespace mrchem {

class H_BM_dia final : public RankTwoTensorOperator<3,3> {
public:
    H_BM_dia(const double *o = 0, const double *k = 0)
            : r_m3(3.0, k),
              r_o(o),
              r_k(k) {
        const double alpha_2 = PHYSCONST::alpha * PHYSCONST::alpha;
        RankZeroTensorOperator &o_x = this->r_o[0];
        RankZeroTensorOperator &o_y = this->r_o[1];
        RankZeroTensorOperator &o_z = this->r_o[2];
        RankZeroTensorOperator &k_x = this->r_k[0];
        RankZeroTensorOperator &k_y = this->r_k[1];
        RankZeroTensorOperator &k_z = this->r_k[2];

        RankTwoTensorOperator<3,3> &h = (*this);
        h[0][0] = -(alpha_2/2.0)*r_m3*(o_y*k_y + o_z*k_z);
        h[0][1] =  (alpha_2/2.0)*r_m3*(o_x*k_y);
        h[0][2] =  (alpha_2/2.0)*r_m3*(o_x*k_z);
        h[1][0] =  (alpha_2/2.0)*r_m3*(o_y*k_x);
        h[1][1] = -(alpha_2/2.0)*r_m3*(o_x*k_x + o_z*k_z);
        h[1][2] =  (alpha_2/2.0)*r_m3*(o_y*k_z);
        h[2][0] =  (alpha_2/2.0)*r_m3*(o_z*k_x);
        h[2][1] =  (alpha_2/2.0)*r_m3*(o_z*k_y);
        h[2][2] = -(alpha_2/2.0)*r_m3*(o_x*k_x + o_y*k_y);
    }
    ~H_BM_dia() { }

protected:
    PowerPotential r_m3;
    PositionOperator r_o;
    PositionOperator r_k;
};

} //namespace mrchem
