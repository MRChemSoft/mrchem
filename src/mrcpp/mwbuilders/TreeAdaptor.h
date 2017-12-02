#pragma once

#include "mrcpp_declarations.h"
#include "TelePrompter.h"

template<int D>
class TreeAdaptor {
public:
    TreeAdaptor(int ms) : maxScale(ms) { }
    virtual ~TreeAdaptor() { }

    void setMaxScale(int ms) { this->maxScale = ms; }

    void splitNodeVector(MWNodeVector &out, MWNodeVector &inp) const {
        for (int n = 0; n < inp.size(); n++) {
            MWNode<D> &node = *inp[n];
            // Can be BranchNode in operator application
            if (node.isBranchNode()) continue;
            if (node.getScale() + 2 > this->maxScale) continue;
            if (splitNode(node)) {
                node.createChildren();
                for (int i = 0; i < node.getNChildren(); i++) {
                    out.push_back(&node.getMWChild(i));
                }
            }
        }
    }

protected:
    int maxScale;

    virtual bool splitNode(const MWNode<D> &node) const { return false; }
};

