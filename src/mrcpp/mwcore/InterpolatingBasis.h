#pragma once

#include "ScalingBasis.h"

class InterpolatingBasis : public ScalingBasis {
public:
    InterpolatingBasis(int k) : ScalingBasis(k, Interpol) {
        initScalingBasis();
        calcQuadratureValues();
        calcCVMaps();
    }
    virtual ~InterpolatingBasis() { }

    void initScalingBasis();
    void calcQuadratureValues();
    void calcCVMaps();
};

