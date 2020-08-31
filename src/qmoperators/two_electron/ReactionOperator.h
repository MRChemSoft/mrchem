#pragma once

#include "ReactionPotential.h"
#include "qmoperators/RankZeroTensorOperator.h"

/** @class ReactionOperator
 *
 * @brief Operator containing a single ReactionPotential
 *
 * This class is a simple TensorOperator realization of @class ReactionPotential.
 *
 */

namespace mrchem {
class SCRF;

class ReactionOperator final : public RankZeroTensorOperator {
public:
    ReactionOperator(OrbitalVector_p Phi_p, SCRF help, bool var = false) {
        potential = std::make_shared<ReactionPotential>(Phi_p, help, var);
        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
    }

    ComplexDouble trace(OrbitalVector &Phi) { return RankZeroTensorOperator::trace(Phi); }

    double getTotalEnergy() { return this->potential->getTotalEnergy(); }
    double getElectronicEnergy() { return this->potential->getElectronicEnergy(); }
    double getNuclearEnergy() { return this->potential->getNuclearEnergy(); }
    bool getRunVariational() { return this->potential->getRunVariational(); }
    SCRF getHelper() { return this->potential->getHelper(); }
    std::shared_ptr<ReactionPotential> getPotential() { return this->potential; }
    void updateMOResidual(double const err_t) { this->potential->updateMOResidual(err_t); }

    QMFunction &getCurrentReactionPotential() { return this->potential->getCurrentReactionPotential(); }
    QMFunction &getPreviousReactionPotential() { return this->potential->getPreviousReactionPotential(); }
    QMFunction &getCurrentDifferenceReactionPotential() {
        return this->potential->getCurrentDifferenceReactionPotential();
    }

    QMFunction &getCurrentGamma() { return this->potential->getCurrentGamma(); }
    QMFunction &getPreviousGamma() { return this->potential->getPreviousGamma(); }
    QMFunction &getCurrentDifferenceGamma() { return this->potential->getCurrentDifferenceGamma(); }

private:
    std::shared_ptr<ReactionPotential> potential{nullptr};
};

} // namespace mrchem
