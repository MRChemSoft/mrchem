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
    ReactionOperator(mrcpp::PoissonOperator *P, mrcpp::DerivativeOperator<3> *D, Cavity *C, const Nuclei &nuc, OrbitalVector *Phi)
            : potential(new ReactionPotential(P, D, C, nuc, Phi)) {
        RankZeroTensorOperator &J = (*this);
        J = *this->potential;
    }
    
    ~ReactionOperator() {
        if (this->potential != nullptr) delete this->potential;
    }


    ComplexDouble trace(OrbitalVector &Phi) { return 0.5 * RankZeroTensorOperator::trace(Phi); }
    
    double &get_tot_Energy(){return this->potential->get_tot_Energy();}
    double &get_e_Energy(){return this->potential->get_e_Energy();}
    double &get_nuc_Energy(){return this->potential->get_nuc_Energy();}
   
private:
    ReactionPotential *potential;
};

} //namespace mrchem
