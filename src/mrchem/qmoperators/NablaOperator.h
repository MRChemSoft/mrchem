#pragma once

#include "QMDerivative.h"
#include "QMTensorOperator.h"

class QMNabla : public QMDerivative {
public:
    QMNabla(int d, DerivativeOperator<3> &D) : QMDerivative(d, D) { }
    virtual ~QMNabla() { }

    virtual void setup(double prec) { setApplyPrec(prec); }
    virtual void clear() { clearApplyPrec(); }

protected:
    virtual bool isReal() const { return true; }
    virtual bool isImag() const { return false; }
};

class NablaOperator : public RankOneTensorOperator<3> {
public:
    NablaOperator(DerivativeOperator<3> &D) : d_x(0, D), d_y(1, D), d_z(2, D) {
        initializeTensorOperator();
    }
    virtual ~NablaOperator() { }

protected:
    QMNabla d_x;
    QMNabla d_y;
    QMNabla d_z;

    void initializeTensorOperator() {
        RankOneTensorOperator<3> &h = (*this);
        h[0] = d_x;
        h[1] = d_y;
        h[2] = d_z;
    }
};


