#pragma once

#include "MWTree.h"
#include "MWNode.h"
#include "TreeIterator.h"

template<int D>
class LebesgueIterator: public TreeIterator<D> {
public:
    LebesgueIterator(MWTree<D> *tree, int dir = TopDown):
        TreeIterator<D>(dir) {
        this->init(tree);
    }
    virtual ~LebesgueIterator() {}
protected:
    int getChildIndex(int i) const { return i; }
};

