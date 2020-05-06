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
    ReactionOperator(std::shared_ptr<mrcpp::PoissonOperator> P,
                     std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                     int history) {
        potential = std::make_shared<ReactionPotential>(P, D, history);
        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
    }

    ComplexDouble trace(OrbitalVector &Phi) { return RankZeroTensorOperator::trace(Phi); }

    double getTotalEnergy() { return this->potential->getTotalEnergy(); }
    double getElectronicEnergy() { return this->potential->getElectronicEnergy(); }
    double getNuclearEnergy() { return this->potential->getNuclearEnergy(); }
    bool getRunVariational() { return this->potential->getRunVariational(); }
    void setRunVariational(bool var) { this->potential->setRunVariational(var); }
    std::shared_ptr<SCRF> getHelper() { return this->potential->getHelper(); }
    void setHelper(std::shared_ptr<SCRF> helper) { this->potential->setHelper(helper); }
    std::shared_ptr<ReactionPotential> getPotential() { return this->potential; }
    void updateTotalDensity(OrbitalVector Phi,double prec) { this->potential->updateTotalDensity(Phi, prec); }

private:
    std::shared_ptr<ReactionPotential> potential{nullptr};
};

} // namespace mrchem
