#include "RankTwoTensorOperator.h"
#include "Orbital.h"
#include "orbital_utils.h"

namespace mrchem {

template<int I, int J>
ComplexMatrix RankTwoTensorOperator<I,J>::operator()(Orbital bra, Orbital ket) {
    RankTwoTensorOperator<I,J> &O = *this;
    ComplexMatrix out(I,J);
    for (int i = 0; i < I; i++) {
        out.row(i) = O[i](bra, ket);
    }
    return out;
}

template<int I, int J>
ComplexMatrix RankTwoTensorOperator<I,J>::trace(OrbitalVector &phi) {
    RankTwoTensorOperator<I,J> &O = *this;
    ComplexMatrix out(I,J);
    for (int i = 0; i < I; i++) {
        out.row(i) = O[i].trace(phi);
    }
    return out;
}

template<int I, int J>
ComplexMatrix RankTwoTensorOperator<I,J>::trace(OrbitalVector &phi, OrbitalVector &x, OrbitalVector &y) {
    RankTwoTensorOperator<I,J> &O = *this;
    ComplexMatrix out(I,J);
    for (int i = 0; i < I; i++) {
        out.row(i) = O[i].trace(phi, x, y);
    }
    return out;
}

} //namespace mrchem

template class mrchem::RankTwoTensorOperator<3,3>;
