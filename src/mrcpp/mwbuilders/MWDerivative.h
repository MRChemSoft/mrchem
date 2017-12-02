#pragma once

#include "TreeBuilder.h"
#include "CopyAdaptor.h"
#include "DefaultCalculator.h"
#include "DerivativeCalculator.h"
#include "DerivativeOperator.h"
#include "Timer.h"

template<int D>
class MWDerivative {
public:
    MWDerivative(int ms = MaxScale) : maxScale(ms) { }
    virtual ~MWDerivative() { }

    int getMaxScale() const { return this->maxScale; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    void operator()(FunctionTree<D> &out,
                    DerivativeOperator<D> &oper,
                    FunctionTree<D> &inp,
                    int dir = -1) const {
        TreeBuilder<D> builder;

        int bw[D]; // Operator bandwidth in [x,y,z]
        for (int d = 0; d < D; d++) bw[d] = 0;

        // Copy input tree plus bandwidth in operator direction
        Timer pre_t;
        oper.calcBandWidths(1.0); // Fixed 0 or 1 for derivatives
        bw[dir] = oper.getMaxBandWidth();
        CopyAdaptor<D> pre_adaptor(inp, this->maxScale, bw);
        DefaultCalculator<D> pre_calculator;
        builder.build(out, pre_calculator, pre_adaptor, -1);
        pre_t.stop();

        // Apply operator on fixed expanded grid
        TreeAdaptor<D> apply_adaptor(this->maxScale);
        DerivativeCalculator<D> apply_calculator(dir, oper, inp);
        builder.build(out, apply_calculator, apply_adaptor, 0);

        Timer post_t;
        oper.clearBandWidths();
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        inp.deleteGenerated();
        post_t.stop();

        println(10, "Time pre operator   " << pre_t);
        println(10, "Time post operator  " << post_t);
        println(10, std::endl);
    }
protected:
    int maxScale;
};

