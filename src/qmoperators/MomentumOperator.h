#pragma once

#include "QMOperator.h"
#include "RankOneTensorOperator.h"

namespace mrchem {

class QMMomentum : public QMOperator {
public:
    QMMomentum(int d, mrcpp::DerivativeOperator<3> &D);
    virtual ~QMMomentum() { }

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

class MomentumOperator : public RankOneTensorOperator<3> {
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
    virtual ~MomentumOperator() { }

protected:
    QMMomentum p_x;
    QMMomentum p_y;
    QMMomentum p_z;
};

} //namespace mrchem
