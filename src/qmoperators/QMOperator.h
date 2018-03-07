#pragma once

#include "qmfunctions.h"

/** \brief Quantum mechanical operators

    Base class to handle operators and their application in the QM sense

  */
namespace mrchem {

class QMOperator {
public:
    QMOperator() : apply_prec(-1.0) { }
    QMOperator(const QMOperator &oper) : apply_prec(oper.apply_prec) { }
    QMOperator& operator=(const QMOperator &inp) { this->apply_prec = inp.apply_prec; return *this; }
    virtual ~QMOperator() { }

    void setApplyPrec(double prec);
    void clearApplyPrec() { this->apply_prec = -1.0; }
    double getApplyPrec() const { return this->apply_prec; }

    virtual void setup(double prec) = 0;
    virtual void clear() = 0;

    virtual Orbital operator()(Orbital inp) = 0;
    virtual Orbital dagger(Orbital inp) = 0;

    virtual OrbitalVector operator()(OrbitalVector inp);
    virtual OrbitalVector dagger(OrbitalVector inp);

    virtual ComplexDouble operator()(Orbital bra, Orbital ket);
    virtual ComplexDouble dagger(Orbital bra, Orbital ket);

    virtual ComplexMatrix operator()(OrbitalVector bra, OrbitalVector ket);
    virtual ComplexMatrix dagger(OrbitalVector bra, OrbitalVector ket);

    virtual ComplexDouble trace(OrbitalVector inp);

protected:
    double apply_prec;
};

} //namespace mrchem
