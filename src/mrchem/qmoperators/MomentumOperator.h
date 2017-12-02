#pragma once

#include "QMDerivative.h"
#include "QMTensorOperator.h"

class QMMomentum : public QMDerivative {
public:
    QMMomentum(int d, DerivativeOperator<3> &D) : QMDerivative(d, D) { }
    virtual ~QMMomentum() { }

    virtual void setup(double prec) { setApplyPrec(prec); }
    virtual void clear() { clearApplyPrec(); }

protected:
    virtual bool isReal() const { return false; }
    virtual bool isImag() const { return true; }
};

class MomentumOperator : public RankOneTensorOperator<3> {
public:
    MomentumOperator(DerivativeOperator<3> &D) : p_x(0, D), p_y(1, D), p_z(2, D) {
        initializeTensorOperator();
    }
    virtual ~MomentumOperator() { }

protected:
    QMMomentum p_x;
    QMMomentum p_y;
    QMMomentum p_z;

    void initializeTensorOperator() {
        RankOneTensorOperator<3> &h = (*this);
        h[0] = p_x;
        h[1] = p_y;
        h[2] = p_z;
    }
};


