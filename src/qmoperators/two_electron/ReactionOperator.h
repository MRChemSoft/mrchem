#pragma once

#include "ReactionPotential.h"
//#include "ReactionPotentialD1.h"
//#include "ReactionPotentialD2.h"
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
                     double eps_i = 1.0,
                     double eps_o = 2.0,
                     bool islin = false) {
        potential = std::make_shared<ReactionPotential>(P, D, C, nuc, Phi, history, eps_i, eps_o, islin);
        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &J = (*this);
        J = potential;
    }

    ~ReactionOperator() override = default;

    ComplexDouble trace(OrbitalVector &Phi) { return RankZeroTensorOperator::trace(Phi); }

    double &getTotalEnergy() { return this->potential->getTotalEnergy(); }
    double &getElectronicEnergy() { return this->potential->getElectronicEnergy(); }
    double &getNuclearEnergy() { return this->potential->getNuclearEnergy(); }
    QMFunction &getGamma() { return this->potential->getGamma(); }
    QMFunction &getGammanp1() { return this->potential->getGammanp1(); }
    void setGamma(QMFunction new_gamma) { this->potential->setGamma(new_gamma); }
    void setGammanp1(QMFunction new_gamma) { this->potential->setGammanp1(new_gamma); }

private:
    std::shared_ptr<ReactionPotential> potential{nullptr};
};

} // namespace mrchem