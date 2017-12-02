#pragma once

#pragma GCC system_header
#include <Eigen/Dense>

#include "ScalingBasis.h"

class LegendreBasis : public ScalingBasis {
public:
    LegendreBasis(int k) : ScalingBasis(k, Legendre) {
        initScalingBasis();
        calcQuadratureValues();
        calcCVMaps();
    }
    virtual ~LegendreBasis() { }

    void initScalingBasis();
    void calcQuadratureValues();
    void calcCVMaps();
};

