#pragma once

#include "qmfunctions.h"
#include "qmoperators.h"

/** \brief Quantum mechanical operators

    Base class to handle operators and their application in the QM sense

  */
namespace mrchem {

class QMOperator {
public:
    QMOperator() : apply_prec(-1.0) { }
    virtual ~QMOperator() { }

    double prec() { return this->apply_prec; }

    friend RankZeroTensorOperator;

protected:
    double apply_prec;

    void setApplyPrec(double prec);
    void clearApplyPrec() { this->apply_prec = -1.0; }

    bool isSetup(double prec);

    virtual void setup(double prec) = 0;
    virtual void clear() = 0;

    virtual Orbital apply(Orbital inp) = 0;
    virtual Orbital dagger(Orbital inp) = 0;
};

} //namespace mrchem
