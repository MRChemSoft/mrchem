#pragma once

#include "GreensKernel.h"

class HelmholtzKernel: public GreensKernel {
public:
    HelmholtzKernel(double m, double eps, double r_min, double r_max)
            : GreensKernel(eps, r_min, r_max),
              mu(m) {
        initializeKernel();
    }
    virtual ~HelmholtzKernel() { }
protected:
    const double mu; /**< exponent */
    virtual void initializeKernel();
};

