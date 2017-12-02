#pragma once

#include "MWNode.h"
#include "mrcpp_declarations.h"

template<int D>
class TreeCalculator {
public:
    TreeCalculator() { }
    virtual ~TreeCalculator() { }

    virtual MWNodeVector* getInitialWorkVector(MWTree<D> &tree) const {
        return tree.copyEndNodeTable();
    }

    virtual void calcNodeVector(MWNodeVector &nodeVec) {
#pragma omp parallel shared(nodeVec)
{
        int nNodes = nodeVec.size();
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D> &node = *nodeVec[n];
            calcNode(node);
        }
}
        postProcess();
    }
protected:
    virtual void calcNode(MWNode<D> &node) = 0;
    virtual void postProcess() { }
};

