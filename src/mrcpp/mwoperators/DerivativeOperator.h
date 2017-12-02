#pragma once

#include "MWOperator.h"

template<int D>
class DerivativeOperator : public MWOperator {
public:
    DerivativeOperator(const MultiResolutionAnalysis<D> &mra) 
            : MWOperator(mra.getOperatorMRA()) { }
    virtual ~DerivativeOperator() { }
};

