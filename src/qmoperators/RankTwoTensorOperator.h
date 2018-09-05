#pragma once

#include "TensorOperator.h"
#include "RankOneTensorOperator.h"

namespace mrchem {

/** @class RankTwoTensorOperator
 *
 *  @brief Matrix of RankZeroTensorOperator
 *
 * This class provides a base for all matrix operators (implemented as a vector
 * of vectors), and implements some simple collective operations returning matrix
 * quantities.
 *
 */

template<int I, int J>
class RankTwoTensorOperator : public TensorOperator<I, RankOneTensorOperator<J> > {
public:
    ComplexMatrix operator()(Orbital bra, Orbital ket);
    ComplexMatrix trace(OrbitalVector &phi);
    ComplexMatrix trace(OrbitalVector &phi, OrbitalVector &x, OrbitalVector &y);
};

} //namespace mrchem
