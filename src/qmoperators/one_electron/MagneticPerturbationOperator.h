#pragma once
#include "tensor/RankZeroOperator.h"
#include "H_B_dip.h"
#include "H_BB_dia.h"

namespace mrchem {

/** @brief Magnetic perturbation operator for uniform field B = (0, 0, B3)
 *
 * W = -(B3/c)*H_B_dip[2] - (B3^2/(2c^2))*H_BB_dia[2][2]
 *   = -(B3/2c)*Lz - (B3^2/8c^2)*(x^2+y^2)
 *
 * Arises from expanding -1/2*nabla_B^2 = -1/2*nabla^2 - W,
 * i.e. T_B = T_0 + W.  The SCF resolvent uses energies shifted by W. */
class MagneticPerturbationOperator final : public RankZeroOperator {
public:
    MagneticPerturbationOperator(double B3, double c,
                                  std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                                  const mrcpp::Coord<3> &o)
            : MagneticPerturbationOperator(B3, c, H_B_dip(D, o), H_BB_dia(o)) {}

private:
    MagneticPerturbationOperator(double B3, double c, H_B_dip para, H_BB_dia dia) {
        RankZeroOperator &W = (*this);
        W = -(B3 / c) * para[2]
          - (B3 * B3 / (2.0 * c * c)) * dia[2][2];
        W.name() = "W_B";
    }
};

} // namespace mrchem