#pragma once
#include "tensor/RankZeroOperator.h"
#include "H_BB_dia.h"
#include <MRCPP/operators/DerivativeOperator.h>

namespace mrchem {

/**
 * Full magnetic perturbation operator for a uniform field B = (0, 0, B3)
 * in symmetric gauge (origin o):
 *
 *   W = (i·B3)/(2c) · Lz  –  (B3²)/(8c²) · (x² + y²)
 *
 * where  Lz = x·∂y – y·∂x.
 *
 * Para part is imaginary; dia part is real.  operator() assembles the
 * complex result directly, without relying on RankZeroOperator's
 * linear_combination machinery (which does not honour the imaginary flag).
 *
 * Input orbitals are expected to be complex (call orbital::CopyToComplex
 * before the SCF loop).
 */
class MagneticPerturbation final : public RankZeroOperator {
public:
    MagneticPerturbation(double B3, double c,
                         std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                         const mrcpp::Coord<3> &o);

    void   setup(double prec);
    void   clear();

    Orbital       operator()(Orbital inp);
    OrbitalVector operator()(OrbitalVector &inp);

private:
    double para_coeff_;                              ///< B3/(2c)
    mrcpp::Coord<3> origin_;
    std::shared_ptr<mrcpp::DerivativeOperator<3>> D_;
    RankZeroOperator dia_part_;                      ///< –(B3²/2c²)·H_BB_dia[2][2]
    double prec_{-1.0};
};

} // namespace mrchem