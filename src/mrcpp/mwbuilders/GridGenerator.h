#pragma once

#include "TreeBuilder.h"
#include "DefaultCalculator.h"
#include "AnalyticAdaptor.h"
#include "CopyAdaptor.h"

template<int D>
class GridGenerator {
public:
    GridGenerator(int ms = MaxScale)
        : maxScale(ms) { }
    virtual ~GridGenerator() { }

    int getMaxScale() const { return this->maxScale; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    void operator()(FunctionTree<D> &out,
                    const RepresentableFunction<D> &inp,
                    int maxIter = -1) const {
        TreeBuilder<D> builder;
        AnalyticAdaptor<D> adaptor(inp, this->maxScale);
        DefaultCalculator<D> calculator;
        builder.build(out, calculator, adaptor, maxIter);
        println(10, std::endl);
    }

    void operator()(FunctionTree<D> &out,
                    FunctionTree<D> &inp,
                    int max_iter = -1) const {
        TreeBuilder<D> builder;
        CopyAdaptor<D> adaptor(inp, this->maxScale, 0);
        DefaultCalculator<D> calculator;
        builder.build(out, calculator, adaptor, max_iter);
        println(10, std::endl);
    }

    void operator()(FunctionTree<D> &out,
                    FunctionTreeVector<D> &inp,
                    int max_iter = -1) const {
        TreeBuilder<D> builder;
        CopyAdaptor<D> adaptor(inp, this->maxScale, 0);
        DefaultCalculator<D> calculator;
        builder.build(out, calculator, adaptor, max_iter);
        println(10, std::endl);
    }
protected:
    int maxScale;
};

