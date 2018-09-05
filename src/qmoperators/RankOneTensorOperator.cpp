#include "RankOneTensorOperator.h"
#include "Orbital.h"
#include "orbital_utils.h"

namespace mrchem {

template<int I>
OrbitalVector RankOneTensorOperator<I>::operator()(Orbital phi) {
    RankOneTensorOperator<I> &O = *this;
    OrbitalVector out;
    for (int i = 0; i < I; i++) {
        out.push_back(O[i](phi));
    }
    return out;
}

template<int I>
ComplexVector RankOneTensorOperator<I>::operator()(Orbital bra, Orbital ket) {
    RankOneTensorOperator<I> &O = *this;
    ComplexVector out(I);
    for (int i = 0; i < I; i++) {
        out(i) = O[i](bra, ket);
    }
    return out;
}

template<int I>
ComplexVector RankOneTensorOperator<I>::trace(OrbitalVector &phi) {
    RankOneTensorOperator<I> &O = *this;
    ComplexVector out = ComplexVector::Zero(I);
    for (int i = 0; i < I; i++) {
        out(i) = O[i].trace(phi);
    }
    return out;
}

template<int I>
ComplexVector RankOneTensorOperator<I>::trace(OrbitalVector &phi, OrbitalVector &x, OrbitalVector &y) {
    RankOneTensorOperator<I> &O = *this;
    ComplexVector out = ComplexVector::Zero(I);
    for (int i = 0; i < I; i++) {
        out(i) = O[i].trace(phi, x, y);
    }
    return out;
}

} //namespace mrchem

template class mrchem::RankOneTensorOperator<3>;
