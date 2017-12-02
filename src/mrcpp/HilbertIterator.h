#pragma once

#include "MWTree.h"
#include "MWNode.h"
#include "TreeIterator.h"
#include "HilbertPath.h"

template<int D>
class HilbertIterator: public TreeIterator<D> {
public:
    HilbertIterator(MWTree<D> *tree, int dir = TopDown)
            : TreeIterator<D>(dir) {
        this->init(tree);
    }
    virtual ~HilbertIterator() { }

protected:
    int getChildIndex(int i) const {
        const MWNode<D> &node = *this->state->node;
        const HilbertPath<D> &h = node.getHilbertPath();
        return h.getZIndex(i);
    }
};

