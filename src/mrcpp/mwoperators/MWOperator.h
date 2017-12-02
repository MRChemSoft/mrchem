#pragma once

#include <vector>

#include "MultiResolutionAnalysis.h"
#include "OperatorTree.h"

class MWOperator {
public:
    MWOperator(MultiResolutionAnalysis<2> mra) : oper_mra(mra) { }
    virtual ~MWOperator() { this->clear(true); }

    int size() const { return this->oper_exp.size(); }
    void push_back(OperatorTree *oper) { this->oper_exp.push_back(oper); } 
    void clear(bool dealloc = false);

    int getMaxBandWidth(int depth = -1) const;
    const Eigen::VectorXi &getMaxBandWidths() const { return this->band_max; }

    void calcBandWidths(double prec);
    void clearBandWidths();

    OperatorTree &getComponent(int i);
    const OperatorTree &getComponent(int i) const;

    OperatorTree *operator[](int i) { return this->oper_exp[i]; }
    const OperatorTree *operator[](int i) const { return this->oper_exp[i]; }

protected:
    MultiResolutionAnalysis<2> oper_mra;
    OperatorTreeVector oper_exp;
    Eigen::VectorXi band_max;
};

