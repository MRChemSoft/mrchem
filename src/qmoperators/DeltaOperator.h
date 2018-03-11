#pragma once

#include "QMPotential.h"
#include "RankZeroTensorOperator.h"

namespace mrchem {

class QMDelta : public QMPotential {
public:
    QMDelta(const double *o, double expo);
    virtual ~QMDelta() { }

protected:
    mrcpp::GaussFunc<3> func;

    virtual void setup(double prec);
    virtual void clear();
};

class DeltaOperator : public RankZeroTensorOperator {
public:
    DeltaOperator(const double *o = 0, double expo = 1.0e6)
            : delta(0, expo) {
        RankZeroTensorOperator &d = (*this);
        d = delta;
    }
    virtual ~DeltaOperator() { }

protected:
    QMDelta delta;
};

} //namespace mrchem
