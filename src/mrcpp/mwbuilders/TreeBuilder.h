#pragma once

#include "mrcpp_declarations.h"

template<int D>
class TreeBuilder {
public:
    TreeBuilder() { }
    virtual ~TreeBuilder() { }

    void build(MWTree<D> &tree,
               TreeCalculator<D> &calculator,
               TreeAdaptor<D> &adaptor,
               int maxIter) const;

    int clear(MWTree<D> &tree,
              TreeCalculator<D> &calculator,
              TreeAdaptor<D> &adaptor) const;
protected:
    double calcScalingNorm(const MWNodeVector &vec) const;
    double calcWaveletNorm(const MWNodeVector &vec) const;
};

