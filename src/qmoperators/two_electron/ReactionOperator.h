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
    ReactionOperator(mrcpp::PoissonOperator *P,
                     mrcpp::DerivativeOperator<3> *D,
                     Cavity *C,
                     const Nuclei &nuc,
                     OrbitalVector *Phi,
                     int history,
                     double eps_i = 1.0,
                     double eps_o = 2.0,
                     bool islin = false)
            : potential(new ReactionPotential(P, D, C, nuc, Phi, history, eps_i, eps_o, islin)) {
        RankZeroTensorOperator &J = (*this);
        J = *this->potential;
    }

    ~ReactionOperator() {
        if (this->potential != nullptr) delete this->potential;
    }

    ComplexDouble trace(OrbitalVector &Phi) { return RankZeroTensorOperator::trace(Phi); }

    double &getTotalEnergy() { return this->potential->getTotalEnergy(); }
    double &getElectronicEnergy() { return this->potential->getElectronicEnergy(); }
    double &getNuclearEnergy() { return this->potential->getNuclearEnergy(); }
    QMFunction &getGamma() { return this->potential->getGamma(); }
    QMFunction &getGammanp1() { return this->potential->getGammanp1(); }
    void setGamma(QMFunction new_gamma) { this->potential->setGamma(new_gamma); }

private:
    ReactionPotential *potential;
};

} // namespace mrchem
