#pragma once

#include "GreensKernel.h"
#include "GaussFunc.h"

class IdentityKernel : public GreensKernel {
public:
    IdentityKernel(double eps)
            : GreensKernel(eps, -1.0, -1.0) {
        initializeKernel();
    }
    virtual ~IdentityKernel() { }
protected:
    virtual void initializeKernel() {
        double alpha = sqrt(1.0/this->epsilon);
        double coef = pow(alpha/pi, 1.0/2.0);
        GaussFunc<1> gFunc(alpha, coef);
        this->append(gFunc);
    }
};

