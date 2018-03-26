#pragma once

#include "RankZeroTensorOperator.h"
#include "CoulombPotential.h"

namespace mrchem {

class CoulombOperator final : public RankZeroTensorOperator {
public:
    CoulombOperator(mrcpp::PoissonOperator &P, OrbitalVector &Phi)
            : potential(0) {
        this->potential = new CoulombPotential(P, Phi);

        RankZeroTensorOperator &J = (*this);
        J = *potential;
    }
    ~CoulombOperator() { delete this->potential; }

    ComplexDouble trace(OrbitalVector &Phi) { return 0.5*RankZeroTensorOperator::trace(Phi); }

protected:
    QMPotential *potential;
};

} //namespace mrchem
