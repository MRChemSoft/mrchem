#pragma once

#include "QMPotential.h"
#include "RankZeroTensorOperator.h"
#include "NuclearFunction.h"

namespace mrchem {

class DistancePotential final : public QMPotential {
public:
    DistancePotential(double pow, const double *r_K, double S);
    ~DistancePotential() { }

    void setup(double prec);
    void clear();

protected:
    const double power;
    NuclearFunction func;
};

class DistanceOperator : public RankZeroTensorOperator {
public:
    DistanceOperator(double pow, const double *R = 0, double S = 1.0e-7)
            : r_pow(pow, R, S) {
        RankZeroTensorOperator &v = (*this);
        v = r_pow;
    }
    ~DistanceOperator() { }

protected:
    DistancePotential r_pow;
};


} //namespace mrchem
