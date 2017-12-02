#pragma once

#include "TreeAdaptor.h"
#include "FunctionTreeVector.h"
#include "constants.h"

template<int D>
class CopyAdaptor : public TreeAdaptor<D> {
public:
    CopyAdaptor(FunctionTree<D> &t, int ms, int *bw)
            : TreeAdaptor<D>(ms) {
        setBandWidth(bw);
        tree_vec.push_back(&t);
    }
    CopyAdaptor(FunctionTreeVector<D> &t, int ms, int *bw)
            : TreeAdaptor<D>(ms),
              tree_vec(t) {
        setBandWidth(bw);
    }
    virtual ~CopyAdaptor() {
    }

protected:
    int bandWidth[D];
    FunctionTreeVector<D> tree_vec;

    void setBandWidth(int *bw) {
        for (int d = 0; d < D; d++) {
            if (bw != 0) {
                this->bandWidth[d] = bw[d];
            } else {
                this->bandWidth[d] = 0;
            }
        }
    }

    bool splitNode(const MWNode<D> &node) const {
        const NodeIndex<D> &idx = node.getNodeIndex();
        for (int c = 0; c < node.getTDim(); c++) {
            NodeIndex<D> cIdx(idx, c);
            for (int d = 0; d < D; d++) {
                for (int bw = -this->bandWidth[d]; bw <= this->bandWidth[d]; bw++) {
                    NodeIndex<D> bwIdx(cIdx);
                    bwIdx.getTranslation()[d] += bw;
                    for (int i = 0; i < this->tree_vec.size(); i++) {
                        const FunctionTree<D> &func_i = tree_vec.getFunc(i);
                        const MWNode<D> *node_i = func_i.findNode(bwIdx);
                        if (node_i != 0) return true;
                    }
                }
            }
        }
        return false;
    }
};

