#pragma once

#include "NuclearGradientOperator.h"
#include "PositionOperator.h"
#include "qmoperators/RankTwoTensorOperator.h"

namespace mrchem {

class H_MB_dia final : public RankTwoTensorOperator<3, 3> {
public:
    H_MB_dia(const mrcpp::Coord<3> &o, const mrcpp::Coord<3> &k, double c)
            : r_o(o)
            , r_rm3(1.0, k, c) {
        const double alpha_2 = PHYSCONST::alpha * PHYSCONST::alpha;
        RankZeroTensorOperator &o_x = this->r_o[0];
        RankZeroTensorOperator &o_y = this->r_o[1];
        RankZeroTensorOperator &o_z = this->r_o[2];
        RankZeroTensorOperator &k_x = this->r_rm3[0];
        RankZeroTensorOperator &k_y = this->r_rm3[1];
        RankZeroTensorOperator &k_z = this->r_rm3[2];

        // Invoke operator= to assign *this operator
        RankTwoTensorOperator<3, 3> &h = (*this);
        h[0][0] = -(alpha_2 / 2.0) * (o_y * k_y + o_z * k_z);
        h[0][1] = (alpha_2 / 2.0) * (o_y * k_x);
        h[0][2] = (alpha_2 / 2.0) * (o_z * k_x);
        h[1][0] = (alpha_2 / 2.0) * (o_x * k_y);
        h[1][1] = -(alpha_2 / 2.0) * (o_x * k_x + o_z * k_z);
        h[1][2] = (alpha_2 / 2.0) * (o_z * k_y);
        h[2][0] = (alpha_2 / 2.0) * (o_x * k_z);
        h[2][1] = (alpha_2 / 2.0) * (o_y * k_z);
        h[2][2] = -(alpha_2 / 2.0) * (o_x * k_x + o_y * k_y);
        h[0][0].name() = "h_MB_dia[x,x]";
        h[0][1].name() = "h_MB_dia[x,y]";
        h[0][2].name() = "h_MB_dia[x,z]";
        h[1][0].name() = "h_MB_dia[y,x]";
        h[1][1].name() = "h_MB_dia[y,y]";
        h[1][2].name() = "h_MB_dia[y,z]";
        h[2][0].name() = "h_MB_dia[z,x]";
        h[2][1].name() = "h_MB_dia[z,y]";
        h[2][2].name() = "h_MB_dia[z,z]";
    }

private:
    PositionOperator r_o;
    NuclearGradientOperator r_rm3;
};

} // namespace mrchem
