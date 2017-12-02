#pragma once

#include "MWNode.h"
#include "OperatorTree.h"

class OperatorNode : public MWNode<2> {
public:
    OperatorTree &getOperTree() { return static_cast<OperatorTree &>(*this->tree); }
    OperatorNode &getOperParent() { return static_cast<OperatorNode &>(*this->parent); }
    OperatorNode &getOperChild(int i) { return static_cast<OperatorNode &>(*this->children[i]); }

    const OperatorTree &getOperTree() const { return static_cast<const OperatorTree &>(*this->tree); }
    const OperatorNode &getOperParent() const { return static_cast<const OperatorNode &>(*this->parent); }
    const OperatorNode &getOperChild(int i) const { return static_cast<const OperatorNode &>(*this->children[i]); }

    void createChildren() {
        MWNode<2>::createChildren();
        this->clearIsEndNode();
    }
    void genChildren() {
        MWNode<2>::createChildren();
        this->clearIsEndNode();
        this->giveChildrenCoefs();
    }
    void deleteChildren() {
        MWNode<2>::deleteChildren();
        this->setIsEndNode();
    }

    friend class MWNode<2>;
    friend class OperatorTree;
    friend class SerialOperatorTree;

protected:
    OperatorNode() : MWNode<2>() { }
    virtual ~OperatorNode() { }

    double calcComponentNorm(int i) const;

    void dealloc() {
        this->tree->decrementNodeCount(this->getScale());
        this->tree->getSerialTree()->deallocNodes(this->getSerialIx());
    }
};

