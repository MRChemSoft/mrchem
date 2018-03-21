#pragma once

#include "QMPotential.h"
#include "RankZeroTensorOperator.h"

namespace mrchem {

class QMDelta final : public QMPotential {
public:
    QMDelta(const double *o, double expo);
    ~QMDelta() { }

protected:
    mrcpp::GaussFunc<3> func;

    void setup(double prec);
    void clear();
};

class DeltaOperator final : public RankZeroTensorOperator {
public:
    DeltaOperator(const double *o = 0, double expo = 1.0e6)
            : delta(0, expo) {
        RankZeroTensorOperator &d = (*this);
        d = delta;
    }
    ~DeltaOperator() { }

protected:
    QMDelta delta;
};

} //namespace mrchem
