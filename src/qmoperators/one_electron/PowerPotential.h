#pragma once

#include "QMPotential.h"
#include "RankZeroTensorOperator.h"
#include "NuclearFunction.h"

namespace mrchem {

class QMPower final : public QMPotential {
public:
    QMPower(double pow, const double *r_K, double S);
    ~QMPower() { }

    void setup(double prec);
    void clear();

protected:
    const double power;
    NuclearFunction func;
};

class PowerPotential : public RankZeroTensorOperator {
public:
    PowerPotential(double pow, const double *R = 0, double S = 1.0e-7)
            : r_pow(pow, R, S) {
        RankZeroTensorOperator &v = (*this);
        v = r_pow;
    }
    ~PowerPotential() { }

protected:
    QMPower r_pow;
};


} //namespace mrchem
