#pragma once

#include "QMTensorOperator.h"

namespace mrchem {

template<int I, int J>
class RankTwoTensorOperator : public QMTensorOperator<I, RankOneTensorOperator<J> > {
public:
    ComplexMatrix operator()(Orbital bra, Orbital ket) {
        RankTwoTensorOperator &O = *this;
        ComplexMatrix out(I,J);
        for (int i = 0; i < I; i++) {
            out.row(i) = O[i](bra, ket);
        }
        return out;
    }
    ComplexMatrix trace(OrbitalVector &phi) {
        RankTwoTensorOperator &O = *this;
        ComplexMatrix out(I,J);
        for (int i = 0; i < I; i++) {
            out.row(i) = O[i].trace(phi);
        }
        return out;
    }
    ComplexMatrix trace(OrbitalVector &phi, OrbitalVector &x, OrbitalVector &y) {
        RankTwoTensorOperator &O = *this;
        ComplexMatrix out(I,J);
        for (int i = 0; i < I; i++) {
            out.row(i) = O[i].trace(phi, x, y);
        }
        return out;
    }
};

} //namespace mrchem
