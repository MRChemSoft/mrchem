#pragma once

#include "RankZeroTensorOperator.h"
#include "ExchangePotential.h"

namespace mrchem {

class ExchangeOperator final : public RankZeroTensorOperator {
public:
    ExchangeOperator(mrcpp::PoissonOperator &P, OrbitalVector &Phi)
            : exchange(0) {
        this->exchange = new ExchangePotential(P, Phi, true);

        RankZeroTensorOperator &K = (*this);
        K = *exchange;
    }
    ~ExchangeOperator() { delete this->exchange; }

    void setupInternal(double prec) { this->exchange->setupInternal(prec); }
    void rotate(const ComplexMatrix &U) { this->exchange->rotate(U); }

    ComplexDouble trace(OrbitalVector &Phi) { return 0.5*RankZeroTensorOperator::trace(Phi); }

protected:
    ExchangePotential *exchange;
};

} //namespace mrchem
