#pragma once

#include "QMOperator.h"
#include "QMTensorOperator.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

class QMSpin : public QMOperator {
public:
    QMSpin() : QMOperator(MRA->getMaxScale()) { }
    virtual ~QMSpin() { }

    virtual void setup(double prec) { setApplyPrec(prec); }
    virtual void clear() { clearApplyPrec(); }
};


class QMSpinX : public QMSpin {
public:
    QMSpinX() : QMSpin() { }
    virtual ~QMSpinX() { }

    virtual Orbital* operator()(Orbital &phi);
    virtual Orbital* adjoint(Orbital &phi);

    using QMOperator::operator();
    using QMOperator::adjoint;
};

class QMSpinY : public QMSpin {
public:
    QMSpinY() : QMSpin() { }
    virtual ~QMSpinY() { }

    virtual Orbital* operator()(Orbital &phi);
    virtual Orbital* adjoint(Orbital &phi);

    using QMOperator::operator();
    using QMOperator::adjoint;
};

class QMSpinZ : public QMSpin {
public:
    QMSpinZ() : QMSpin() { }
    virtual ~QMSpinZ() { }

    virtual Orbital* operator()(Orbital &phi);
    virtual Orbital* adjoint(Orbital &phi);

    using QMOperator::operator();
    using QMOperator::adjoint;
};

class SpinOperator : public RankOneTensorOperator<3> {
public:
    SpinOperator() : s_x(), s_y(), s_z(){
        initializeTensorOperator();
    }
    virtual ~SpinOperator() { }

protected:
    QMSpinX s_x;
    QMSpinY s_y;
    QMSpinZ s_z;

    void initializeTensorOperator() {
        RankOneTensorOperator<3> &h = (*this);
        h[0] = s_x;
        h[1] = s_y;
        h[2] = s_z;
    }
};

