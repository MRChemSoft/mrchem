#pragma once
#include "tensor/RankZeroOperator.h"
#include "H_B_dip.h"
#include "H_BB_dia.h"

namespace mrchem {

/** Paramagnetic part: W_para = -(B3/c)*H_B_dip[2] = -(B3/2c)*Lz  [imaginary] */
class MagneticPerturbationPara final : public RankZeroOperator {
public:
    MagneticPerturbationPara(double B3, double c,
                              std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                              const mrcpp::Coord<3> &o) {
        RankZeroOperator &h = (*this);
        h = -(B3 / c) * H_B_dip(D, o)[2];
        h.name() = "W_B_para";
    }
};

/** Diamagnetic part: W_dia = -(B3^2/(2c^2))*H_BB_dia[2][2] = -(B3^2/8c^2)*(x^2+y^2)  [real] */
class MagneticPerturbationDia final : public RankZeroOperator {
public:
    MagneticPerturbationDia(double B3, double c,
                             const mrcpp::Coord<3> &o) {
        RankZeroOperator &h = (*this);
        h = -(B3 * B3 / (2.0 * c * c)) * H_BB_dia(o)[2][2];
        h.name() = "W_B_dia";
    }
};

} // namespace mrchem