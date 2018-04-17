#pragma once

#include "RankZeroTensorOperator.h"
#include "ExchangePotential.h"

/** @class ExchangeOperator
 *
 * @brief Operator containing a single ExchangePotential
 *
 * The operator is defined as the Hartree-Fock exchange arising from a particular
 * set of orbitals. The OrbitalVector defining the operator is fixed throughout the
 * operators life time, but the orbitals themselves are allowed to change in between
 * each application. The internal exchange potentials (the operator applied to it's
 * own orbitals) can be precomputed and stored for fast retrieval.
 */

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
