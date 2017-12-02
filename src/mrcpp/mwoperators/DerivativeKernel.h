#pragma once

#include "GreensKernel.h"
#include "GaussFunc.h"
#include "GaussPoly.h"

class DerivativeKernel : public GreensKernel {
public:
    DerivativeKernel(double eps)
            : GreensKernel(eps, -1.0, -1.0) {
        initializeKernel();
    }
    virtual ~DerivativeKernel() { }
protected:
    virtual void initializeKernel() {
        double alpha = 1.0/this->epsilon;
        double coef = pow(alpha/pi, 1.0/2.0);
        GaussFunc<1> g(alpha, coef);
        GaussPoly<1> dg = g.differentiate(0);
        this->append(dg);
    }
};

