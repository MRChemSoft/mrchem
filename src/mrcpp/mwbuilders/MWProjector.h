#pragma once

#include "TreeBuilder.h"
#include "ProjectionCalculator.h"
#include "WaveletAdaptor.h"
#include "AnalyticFunction.h"
#include "Timer.h"

template<int D>
class MWProjector {
public:
    MWProjector(double pr = -1.0, int ms = MaxScale)
        : prec(pr), maxScale(ms) { }
    virtual ~MWProjector() { }

    double getPrecision() const { return this->prec; }
    int getMaxScale() const { return this->maxScale; }

    void setPrecision(double pr) { this->prec = pr; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    void operator()(FunctionTree<D> &out,
                    std::function<double (const double *r)> func,
                    int maxIter = -1) const {
        AnalyticFunction<D> inp(func);
        (*this)(out, inp, maxIter);
    }
    void operator()(FunctionTree<D> &out,
                    RepresentableFunction<D> &inp,
                    int maxIter = -1) const {
        TreeBuilder<D> builder;
        WaveletAdaptor<D> adaptor(this->prec, this->maxScale);
        ProjectionCalculator<D> calculator(inp);

        builder.build(out, calculator, adaptor, maxIter);

        Timer trans_t;
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        trans_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, std::endl);
    }
protected:
    double prec;
    int maxScale;
};

