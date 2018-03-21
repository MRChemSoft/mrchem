#pragma once

#include "RankZeroTensorOperator.h"
#include "MomentumOperator.h"

namespace mrchem {

class KineticOperator final : public RankZeroTensorOperator {
public:
    KineticOperator(mrcpp::DerivativeOperator<3> &D)
            : p(D) {
        RankZeroTensorOperator &p_x = this->p[0];
        RankZeroTensorOperator &p_y = this->p[1];
        RankZeroTensorOperator &p_z = this->p[2];

        RankZeroTensorOperator &t = (*this);
        t = 0.5*(p_x*p_x + p_y*p_y + p_z*p_z);
    }
    ~KineticOperator() { }

    ComplexDouble operator()(Orbital bra, Orbital ket);
    ComplexDouble dagger(Orbital bra, Orbital ket);

    ComplexMatrix operator()(OrbitalVector &bra, OrbitalVector &ket);
    ComplexMatrix dagger(OrbitalVector &bra, OrbitalVector &ket);

    using RankZeroTensorOperator::operator();
    using RankZeroTensorOperator::dagger;

protected:
    MomentumOperator p;
};

} //namespace mrchem
