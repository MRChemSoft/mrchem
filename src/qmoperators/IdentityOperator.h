#pragma once

#include "QMOperator.h"

namespace mrchem {

class IdentityOperator : public QMOperator {
public:
    IdentityOperator() : QMOperator() { }
    virtual ~IdentityOperator() { }

    virtual void setup(double prec) { setApplyPrec(prec); }
    virtual void clear() { clearApplyPrec(); }

    virtual Orbital operator()(Orbital inp);
    virtual Orbital dagger(Orbital inp);

    //virtual OrbitalVector operator()(OrbitalVector inp);
    //virtual OrbitalVector dagger(OrbitalVector inp);

    virtual ComplexDouble operator()(Orbital bra, Orbital ket);
    virtual ComplexDouble dagger(Orbital bra, Orbital ket);

    //virtual ComplexMatrix operator()(OrbitalVector bra, OrbitalVector ket);
    //virtual ComplexMatrix dagger(OrbitalVector bra, OrbitalVector ket);

    using QMOperator::operator();
    using QMOperator::dagger;

protected:
    //ComplexMatrix calcOverlapMatrix(OrbitalVector bra, OrbitalVector ket);
    //ComplexMatrix calcOverlapMatrix_P(OrbitalVector bra, OrbitalVector ket);
    //ComplexMatrix calcOverlapMatrix_P_H(OrbitalVector bra, OrbitalVector ket);
};

} //namespace mrchem
