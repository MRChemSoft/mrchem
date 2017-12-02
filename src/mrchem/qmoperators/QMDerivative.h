#pragma once

#include "QMOperator.h"

template<int D> class DerivativeOperator;

class QMDerivative : public QMOperator {
public:
    QMDerivative(int dir, DerivativeOperator<3> &d);
    virtual ~QMDerivative() { }

    virtual Orbital* operator() (Orbital &phi_p);
    virtual Orbital* adjoint(Orbital &phi_p);

    using QMOperator::operator();
    using QMOperator::adjoint;

protected:
    const int apply_dir;
    DerivativeOperator<3> *derivative;

    virtual bool isReal() const = 0;
    virtual bool isImag() const = 0;

    void calcRealPart(Orbital &dPhi, Orbital &phi);
    void calcImagPart(Orbital &dPhi, Orbital &phi);
};

