#pragma once

#include "RankZeroTensorOperator.h"

/** @class FockOperator
 *
 * @brief Operator containing the standard SCF operators
 *
 * This is a simple collection of operators used in ground state SCF calculations.
 * The operator is separated into kinetic and potential parts, since the MW way of
 * solving the SCF equations is to invert the kinetic part, and apply the potential
 * part as usual.
 */

namespace mrchem {

class SCFEnergy;
class KineticOperator;
class NuclearOperator;
class CoulombOperator;
class ExchangeOperator;
class XCOperator;

class FockOperator final : public RankZeroTensorOperator {
public:
    FockOperator(KineticOperator  *t = 0,
                 NuclearOperator  *v = 0,
                 CoulombOperator  *j = 0,
                 ExchangeOperator *k = 0,
                 XCOperator       *xc = 0);
    ~FockOperator() { }

    RankZeroTensorOperator& kinetic()   { return this->T; }
    RankZeroTensorOperator& potential() { return this->V; }

    KineticOperator           *getKineticOperator()  { return this->kin;  }
    NuclearOperator           *getNuclearOperator()  { return this->nuc;  }
    CoulombOperator           *getCoulombOperator()  { return this->coul; }
    ExchangeOperator          *getExchangeOperator() { return this->ex;   }
    XCOperator                *getXCOperator()       { return this->xc;   }

    void rotate(const ComplexMatrix &U);

    void setup(double prec);
    void clear();
    void addExternalPotential(RankZeroTensorOperator &O);
    
    SCFEnergy trace(OrbitalVector &Phi, const ComplexMatrix &F);

protected:
    RankZeroTensorOperator T;     ///< Total kinetic energy operator
    RankZeroTensorOperator V;     ///< Total potential energy operator
    RankZeroTensorOperator V_ext; ///< Total external potential
    
    KineticOperator *kin;
    NuclearOperator *nuc;
    CoulombOperator *coul;
    ExchangeOperator *ex;
    XCOperator *xc;
};

} //namespace mrchem
