#pragma once

#include "TreeCalculator.h"

template<int D>
class DefaultCalculator :public TreeCalculator<D> {
public:
    DefaultCalculator() { }
    virtual ~DefaultCalculator() { }

    // Reimplementation without OpenMP, the default is faster this way
    virtual void calcNodeVector(MWNodeVector &nodeVec) {
        int nNodes = nodeVec.size();
        for (int n = 0; n < nNodes; n++) {
            calcNode(*nodeVec[n]);
        }
    }
protected:
    virtual void calcNode(MWNode<D> &node) {
        node.clearHasCoefs();
        node.clearNorms();
    }
};

