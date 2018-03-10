#pragma once

#include "QMTensorOperator.h"
#include "RankZeroTensorOperator.h"
#include "Orbital.h"

namespace mrchem {

template<int I>
class RankOneTensorOperator : public QMTensorOperator<I, RankZeroTensorOperator> {
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
