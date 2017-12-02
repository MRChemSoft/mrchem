#pragma once

#include "ConvolutionOperator.h"
#include "HelmholtzKernel.h"

class HelmholtzOperator : public ConvolutionOperator<3> {
public:
    HelmholtzOperator(const MultiResolutionAnalysis<3> &mra,
                      double m, double pr = -1.0)
            : ConvolutionOperator<3>(mra, pr), mu(m) {
        int oldlevel = TelePrompter::setPrintLevel(0);
        double epsilon = this->prec/10.0;
        double r_min = calcMinDistance(mra, epsilon);
        double r_max = calcMaxDistance(mra);
        HelmholtzKernel helmholtz_kernel(this->mu, epsilon, r_min, r_max);
        // Rescale for application in 3D
        helmholtz_kernel.rescale(3);
        initializeOperator(helmholtz_kernel);
        TelePrompter::setPrintLevel(oldlevel);
    }
    virtual ~HelmholtzOperator() { }

    double getMu() const { return this->mu; }
protected:
    const double mu;
};


