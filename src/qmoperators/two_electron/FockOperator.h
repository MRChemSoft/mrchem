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
class ElectricFieldOperator;

class FockOperator final : public RankZeroTensorOperator {
public:
    FockOperator(KineticOperator         *t = nullptr,
                 NuclearOperator         *v = nullptr,
                 CoulombOperator         *j = nullptr,
                 ExchangeOperator        *k = nullptr,
                 XCOperator             *xc = nullptr,
                 ElectricFieldOperator *ext = nullptr);

    KineticOperator        &kinetic()       { return *this->kin; }
    RankZeroTensorOperator &potential()     { return this->V; }
    RankZeroTensorOperator &perturbation()  { return this->H_1; }

    KineticOperator        *getKineticOperator()  { return this->kin;  }
    NuclearOperator        *getNuclearOperator()  { return this->nuc;  }
    CoulombOperator        *getCoulombOperator()  { return this->coul; }
    ExchangeOperator       *getExchangeOperator() { return this->ex;   }
    XCOperator             *getXCOperator()       { return this->xc;   }
    ElectricFieldOperator  *getExtOperator()      { return this->ext;  }
    
    void setKineticOperator (KineticOperator        *t)  { this->kin = t;  }
    void setNuclearOperator (NuclearOperator        *v)  { this->nuc = v;  }
    void setCoulombOperator (CoulombOperator        *j)  { this->coul = j; }
    void setExchangeOperator(ExchangeOperator       *k)  { this->ex = k;   }
    void setXCOperator      (XCOperator             *xc) { this->xc = xc;  }
    void setExtOperator     (ElectricFieldOperator  *ext){ this->ext = ext;}
    
    void rotate(const ComplexMatrix &U);

    void build();
    void setup(double prec);
    void clear();
    
    SCFEnergy trace(OrbitalVector &Phi, const ComplexMatrix &F);

protected:
    RankZeroTensorOperator V;     ///< Total potential energy operator
    RankZeroTensorOperator H_1;   ///< Perturbation operators

    KineticOperator        *kin;
    NuclearOperator        *nuc;
    CoulombOperator        *coul;
    ExchangeOperator       *ex;
    XCOperator             *xc;
    ElectricFieldOperator  *ext;  ///< Total external potential
};

} //namespace mrchem
