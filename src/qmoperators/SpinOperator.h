#pragma once

#include "QMOperator.h"
#include "RankOneTensorOperator.h"

namespace mrchem {

class QMSpin : public QMOperator {
public:
    QMSpin(int d) : D(d) { }
    virtual ~QMSpin() { }

protected:
    const int D;

    virtual void setup(double prec) { setApplyPrec(prec); }
    virtual void clear() { clearApplyPrec(); }

    virtual Orbital apply(Orbital inp);
    virtual Orbital dagger(Orbital inp);

    using QMOperator::apply;    // Necessary in order to pick up base class 
    using QMOperator::dagger;   // definitions for overloaded functions
};

class SpinOperator : public RankOneTensorOperator<3> {
public:
    SpinOperator()
        : s_x(0),
          s_y(1),
          s_z(2) {
        RankOneTensorOperator &s = *this;
        s[0] = s_x;
        s[1] = s_y;
        s[2] = s_z;
    }
    virtual ~SpinOperator() { }

protected:
    QMSpin s_x;
    QMSpin s_y;
    QMSpin s_z;
};

} //namespace mrchem
