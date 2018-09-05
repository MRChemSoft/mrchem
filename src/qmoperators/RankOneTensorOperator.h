#pragma once

#include "TensorOperator.h"
#include "RankZeroTensorOperator.h"

namespace mrchem {

/** @class RankOneTensorOperator
 *
 *  @brief Vector of RankZeroTensorOperator
 *
 * This class provides a base for all vector operators, and implements some simple
 * collective operations returning vector quantities.
 *
 */

template<int I>
class RankOneTensorOperator : public TensorOperator<I, RankZeroTensorOperator> {
public:
    OrbitalVector operator()(Orbital phi);
    ComplexVector operator()(Orbital bra, Orbital ket);
    ComplexVector trace(OrbitalVector &phi);
    ComplexVector trace(OrbitalVector &phi, OrbitalVector &x, OrbitalVector &y);
};

} //namespace mrchem
