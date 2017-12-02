#pragma once

#include "WaveletAdaptor.h"
#include "DefaultCalculator.h"
#include "mrcpp_declarations.h"

template<int D>
class GridCleaner {
public:
    GridCleaner(double pr = -1.0, int ms = MaxScale)
        : prec(pr), maxScale(ms) { }
    virtual ~GridCleaner() { }

    double getPrecision() const { return this->prec; }
    int getMaxScale() const { return this->maxScale; }

    void setPrecision(double pr) { this->prec = pr; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    int operator()(MWTree<D> &out) const {
        TreeBuilder<D> builder;
        DefaultCalculator<D> calculator;
        WaveletAdaptor<D> adaptor(this->prec, this->maxScale);
        return builder.clear(out, calculator, adaptor);
    }

protected:
    double prec;
    int maxScale;
};

