#pragma once

#include "QMOperator.h"
#include "RankZeroTensorOperator.h"

namespace mrchem {

class QMIdentity : public QMOperator {
public:
    QMIdentity() : QMOperator() { }
    virtual ~QMIdentity() { }

protected:
    virtual void setup(double prec) { setApplyPrec(prec); }
    virtual void clear() { clearApplyPrec(); }

    virtual Orbital apply(Orbital inp);
    virtual Orbital dagger(Orbital inp);

    virtual ComplexDouble apply(Orbital bra, Orbital ket);
    virtual ComplexDouble dagger(Orbital bra, Orbital ket);

    virtual ComplexMatrix apply(OrbitalVector &bra, OrbitalVector &ket);
    virtual ComplexMatrix dagger(OrbitalVector &bra, OrbitalVector &ket);
};

class IdentityOperator : public RankZeroTensorOperator {
public:
    IdentityOperator() {
        RankZeroTensorOperator &h = (*this);
        h = I;
    }
    virtual ~IdentityOperator() { }

protected:
    QMIdentity I;
};

} //namespace mrchem
