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

class ReactionOperator final : public RankZeroTensorOperator {
public:
    ReactionOperator(std::shared_ptr<mrcpp::PoissonOperator> P,
                     std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                     std::shared_ptr<mrchem::Cavity> C,
                     const Nuclei &nuc,
                     std::shared_ptr<mrchem::OrbitalVector> Phi,
                     int history,
                     double eps_i,
                     double eps_o,
                     bool islin) {
        potential = std::make_shared<ReactionPotential>(P, D, C, nuc, Phi, history, eps_i, eps_o, islin);
        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
    }

    ComplexDouble trace(OrbitalVector &Phi) { return RankZeroTensorOperator::trace(Phi); }

    double &getTotalEnergy() { return this->potential->getTotalEnergy(); }
    double &getElectronicEnergy() { return this->potential->getElectronicEnergy(); }
    double &getNuclearEnergy() { return this->potential->getNuclearEnergy(); }
    bool &getRunVariational() { return this->potential->getRunVariational(); }
    QMFunction &getGamma() { return this->potential->getGamma(); }
    QMFunction &getDiffFunc() { return this->potential->getDiffFunc(); }
    QMFunction &getPotential() { return *this->potential; }
    void setDiffFunc(QMFunction new_diff_func) { this->potential->setDiffFunc(new_diff_func); }
    void setRunVariational(bool var) { this->potential->setRunVariational(var); }

private:
    std::shared_ptr<ReactionPotential> potential{nullptr};
};

} // namespace mrchem
