#pragma once

#include "ConvolutionOperator.h"
#include "PoissonKernel.h"

class PoissonOperator : public ConvolutionOperator<3> {
public:
    PoissonOperator(const MultiResolutionAnalysis<3> &mra, double pr = -1.0)
            : ConvolutionOperator<3>(mra, pr) {
        int oldlevel = TelePrompter::setPrintLevel(0);
        double epsilon = this->prec/10.0;
        double r_min = calcMinDistance(mra, epsilon);
        double r_max = calcMaxDistance(mra);
        PoissonKernel poisson_kernel(epsilon, r_min, r_max);
        // Rescale for application in 3D
        poisson_kernel.rescale(3);
        initializeOperator(poisson_kernel);
        TelePrompter::setPrintLevel(oldlevel);
    }
    virtual ~PoissonOperator() { }
};

