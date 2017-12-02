#pragma once

#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "MWFDCalculator.h"

template<int D>
class MWFDOperator {
public:
    MWFDOperator(const MultiResolutionAnalysis<D> &mra,
                 int k, int n = -1, double pr = -1.0)
        : prec(pr),
          maxScale(mra.getMaxScale()),
          diff_order(k),
          approx_order(n) {
        if (this->approx_order < 0)
            this->approx_order = mra.getOrder();
        if (this->diff_order > this->approx_order)
            MSG_FATAL("Derivative order must be lower than approximation order");
        if (mra.getOrder() < this->approx_order)
            MSG_FATAL("Approximation order must be lower than polynomial order");
    }
    virtual ~MWFDOperator() { }

    void setPrecision(double pr) { this->prec = pr; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    void operator()(FunctionTree<D> &out,
                    FunctionTree<D> &inp,
                    int maxIter = -1,
                    int dir = -1) const {
        Timer pre_t;
        TreeBuilder<D> builder;
        WaveletAdaptor<D> adaptor(this->prec, this->maxScale);
        MWFDCalculator<D> calculator(dir, this->diff_order, this->approx_order, inp);
        pre_t.stop();

        builder.build(out, calculator, adaptor, maxIter);

        Timer post_t;
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        inp.deleteGenerated();
        post_t.stop();

        println(10, "Time pre operator   " << pre_t);
        println(10, "Time post operator  " << post_t);
        println(10, std::endl);
    }

protected:
    double prec;
    int maxScale;
    int diff_order;
    int approx_order;
};

