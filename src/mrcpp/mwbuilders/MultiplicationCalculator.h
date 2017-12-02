#pragma once

#include "TreeCalculator.h"
#include "FunctionTreeVector.h"

template<int D> class MWMultiplier;

template<int D>
class MultiplicationCalculator : public TreeCalculator<D> {
public:
    friend class MWMultiplier<D>;
protected:
    MultiplicationCalculator(FunctionTreeVector<D> &inp) : prod_vec(&inp) { }
    virtual ~MultiplicationCalculator() { }

    virtual void calcNode(MWNode<D> &node_o) {
        const NodeIndex<D> &idx = node_o.getNodeIndex();
        double *coefs_o = node_o.getCoefs();
        for (int j = 0; j < node_o.getNCoefs(); j++) {
            coefs_o[j] = 1.0;
        }
        for (int i = 0; i < this->prod_vec->size(); i++) {
            double c_i = this->prod_vec->getCoef(i);
            FunctionTree<D> &func_i = this->prod_vec->getFunc(i);
            // This generates missing nodes
            MWNode<D> node_i = func_i.getNode(idx); // Copy node
            node_i.mwTransform(Reconstruction);
            node_i.cvTransform(Forward);
            const double *coefs_i = node_i.getCoefs();
            for (int j = 0; j < node_i.getNCoefs(); j++) {
                coefs_o[j] *= c_i * coefs_i[j];
            }
        }
        node_o.cvTransform(Backward);
        node_o.mwTransform(Compression);
        node_o.setHasCoefs();
        node_o.calcNorms();
    }

private:
    FunctionTreeVector<D> *prod_vec;
};

