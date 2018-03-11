#pragma once

#include "QMOperator.h"
#include "RankOneTensorOperator.h"

namespace mrchem {

class QMNabla : public QMOperator {
public:
    QMNabla(int d, mrcpp::DerivativeOperator<3> &D);
    virtual ~QMNabla() { }

protected:
    const int apply_dir;
    mrcpp::DerivativeOperator<3> *derivative;

    virtual void setup(double prec) { setApplyPrec(prec); }
    virtual void clear() { clearApplyPrec(); }

    virtual Orbital apply(Orbital inp);
    virtual Orbital dagger(Orbital inp);

    using QMOperator::apply;    // Necessary in order to pick up base class 
    using QMOperator::dagger;   // definitions for overloaded functions
};

class NablaOperator : public RankOneTensorOperator<3> {
public:
    NablaOperator(mrcpp::DerivativeOperator<3> &D)
            : d_x(0, D),
              d_y(1, D),
              d_z(2, D) {
        RankOneTensorOperator<3> &d = (*this);
        d[0] = d_x;
        d[1] = d_y;
        d[2] = d_z;
    }
    virtual ~NablaOperator() { }

protected:
    QMNabla d_x;
    QMNabla d_y;
    QMNabla d_z;
};

} //namespace mrchem
