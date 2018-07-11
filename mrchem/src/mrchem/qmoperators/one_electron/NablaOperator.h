#pragma once

#include "qmoperators/QMOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class QMNabla final : public QMOperator {
public:
    QMNabla(int d, mrcpp::DerivativeOperator<3> &D);
    ~QMNabla() { }

protected:
    const int apply_dir;
    mrcpp::DerivativeOperator<3> *derivative;

    void setup(double prec) { setApplyPrec(prec); }
    void clear() { clearApplyPrec(); }

    Orbital apply(Orbital inp);
    Orbital dagger(Orbital inp);
};

class NablaOperator final : public RankOneTensorOperator<3> {
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
    ~NablaOperator() { }

protected:
    QMNabla d_x;
    QMNabla d_y;
    QMNabla d_z;
};

} //namespace mrchem
