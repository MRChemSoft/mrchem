#pragma once

#include "MWOperator.h"
#include "FunctionTreeVector.h"

template<int D>
class ConvolutionOperator : public MWOperator {
public:
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra, double pr);
    virtual ~ConvolutionOperator();

protected:
    MultiResolutionAnalysis<1> kern_mra;
    FunctionTreeVector<1> kern_exp;
    double prec;

    void initializeOperator(GreensKernel &greens_kernel);
    void clearKernel();

    double calcMinDistance(const MultiResolutionAnalysis<D> &MRA, double epsilon) const;
    double calcMaxDistance(const MultiResolutionAnalysis<D> &MRA) const;
};

