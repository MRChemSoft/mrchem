#pragma once

#include "TreeCalculator.h"
#include "OperatorStatistics.h"

template<int D>
class DerivativeCalculator : public TreeCalculator<D> {
public:
    DerivativeCalculator(int dir, DerivativeOperator<D> &o, FunctionTree<D> &f);
    virtual ~DerivativeCalculator();

    virtual MWNodeVector* getInitialWorkVector(MWTree<D> &tree) const;

protected:
    int applyDir;
    FunctionTree<D> *fTree;
    DerivativeOperator<D> *oper;

    std::vector<Timer> band_t;
    std::vector<Timer> calc_t;
    std::vector<Timer> norm_t;
    OperatorStatistics<D> operStat;

    MWNodeVector makeOperBand(const MWNode<D> &gNode);

    void initTimers();
    void clearTimers();
    void printTimers() const;

    virtual void calcNode(MWNode<D> &node);
    virtual void postProcess() {
        printTimers();
        clearTimers();
        initTimers();
    }

    void applyOperator(OperatorState<D> &os);
    void tensorApplyOperComp(OperatorState<D> &os);
};


