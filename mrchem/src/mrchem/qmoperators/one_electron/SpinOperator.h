#pragma once

#include "qmoperators/QMOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class QMSpin final : public QMOperator {
public:
    QMSpin(int d) : D(d) { }
    ~QMSpin() { }

protected:
    const int D;

    void setup(double prec) { setApplyPrec(prec); }
    void clear() { clearApplyPrec(); }

    Orbital apply(Orbital inp);
    Orbital dagger(Orbital inp);
};

class SpinOperator final : public RankOneTensorOperator<3> {
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
    ~SpinOperator() { }

protected:
    QMSpin s_x;
    QMSpin s_y;
    QMSpin s_z;
};

} //namespace mrchem
