#pragma once

#include "TensorOperator.h"
#include "RankZeroTensorOperator.h"
#include "Orbital.h"

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
    OrbitalVector operator()(Orbital phi) {
        RankOneTensorOperator &O = *this;
        OrbitalVector out;
        for (int i = 0; i < I; i++) {
            out.push_back(O[i](phi));
        }
        return out;
    }
    ComplexVector operator()(Orbital bra, Orbital ket) {
        RankOneTensorOperator &O = *this;
        ComplexVector out(I);
        for (int i = 0; i < I; i++) {
            out(i) = O[i](bra, ket);
        }
        return out;
    }
    ComplexVector trace(OrbitalVector &phi) {
        RankOneTensorOperator &O = *this;
        ComplexVector out = ComplexVector::Zero(I);
        for (int i = 0; i < I; i++) {
            out(i) = O[i].trace(phi);
        }
        return out;
    }
    ComplexVector trace(OrbitalVector &phi, OrbitalVector &x, OrbitalVector &y) {
        RankOneTensorOperator &O = *this;
        ComplexVector out = ComplexVector::Zero(I);
        for (int i = 0; i < I; i++) {
            out(i) = O[i].trace(phi, x, y);
        }
        return out;
    }
};

} //namespace mrchem
