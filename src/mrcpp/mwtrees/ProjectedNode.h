#pragma once

#include "FunctionNode.h"

template<int D> class SerialFunctionTree;

template<int D>
class ProjectedNode: public FunctionNode<D> {
public:
    void createChildren() {
        MWNode<D>::createChildren();
        this->clearIsEndNode();
    }
    void genChildren() {
        if (this->isBranchNode()) MSG_FATAL("Node already has children");
        this->tree->getSerialTree()->allocGenChildren(*this);
        this->setIsBranchNode();
    }
    void deleteChildren() {
        MWNode<D>::deleteChildren();
        this->setIsEndNode();
    }
    friend class SerialFunctionTree<D>;

protected:
    ProjectedNode() : FunctionNode<D>() { }
    virtual ~ProjectedNode() { assert(this->tree == 0); }

    void dealloc() {
        this->tree->decrementNodeCount(this->getScale());
        this->tree->getSerialTree()->deallocNodes(this->getSerialIx());
    }

    void reCompress();
};

