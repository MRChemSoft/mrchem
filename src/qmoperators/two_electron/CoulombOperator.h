#pragma once

#include "RankZeroTensorOperator.h"
#include "CoulombPotential.h"

/** @class CoulombOperator
 *
 * @brief Operator containing a single QMPotential
 *
 * The QMPotential is defined as the Coulomb potential arising from a particular
 * set of orbitals. The OrbitalVector defining the operator is fixed throughout the
 * operators life time, but the orbitals themselves are allowed to change in between
 * each application. When the operator has been setup, it will be fixed until it's
 * cleared and setup again (even if the orbitals change).
 */

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
