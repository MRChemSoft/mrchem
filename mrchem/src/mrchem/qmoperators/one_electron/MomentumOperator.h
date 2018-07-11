#pragma once

#include "qmoperators/QMOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class QMMomentum final : public QMOperator {
public:
    QMMomentum(int d, mrcpp::DerivativeOperator<3> &D);
    ~QMMomentum() { }

protected:
    const int apply_dir;
    mrcpp::DerivativeOperator<3> *derivative;

    void setup(double prec) { setApplyPrec(prec); }
    void clear() { clearApplyPrec(); }

    Orbital apply(Orbital inp);
    Orbital dagger(Orbital inp);
};

class MomentumOperator final : public RankOneTensorOperator<3> {
public:
    MomentumOperator(mrcpp::DerivativeOperator<3> &D)
            : p_x(0, D),
              p_y(1, D),
              p_z(2, D) {
        RankOneTensorOperator<3> &p = (*this);
        p[0] = p_x;
        p[1] = p_y;
        p[2] = p_z;
    }
    ~MomentumOperator() { }

protected:
    QMMomentum p_x;
    QMMomentum p_y;
    QMMomentum p_z;
};

} //namespace mrchem
