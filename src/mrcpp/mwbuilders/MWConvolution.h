#pragma once

#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "ConvolutionCalculator.h"
#include "ConvolutionOperator.h"
#include "Timer.h"

template<int D>
class MWConvolution {
public:
    MWConvolution(double pr = -1.0, int ms = MaxScale, bool ap = false)
        : prec(pr), absPrec(ap), maxScale(ms) { }
    virtual ~MWConvolution() { }

    double getPrecision() const { return this->prec; }
    bool getAbsPrec() const { return this->absPrec; }
    int getMaxScale() const { return this->maxScale; }

    void setPrecision(double pr) { this->prec = pr; }
    void setAbsPrec(bool ap) { this->absPrec = ap; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    void operator()(FunctionTree<D> &out,
                    ConvolutionOperator<D> &oper,
                    FunctionTree<D> &inp,
                    int max_iter = -1) const {
        Timer pre_t;
        oper.calcBandWidths(this->prec);
        WaveletAdaptor<D> adaptor(this->prec, this->maxScale, this->absPrec);
        ConvolutionCalculator<D> calculator(this->prec, oper, inp);
        pre_t.stop();

        TreeBuilder<D> builder;
        builder.build(out, calculator, adaptor, max_iter);

        Timer post_t;
        oper.clearBandWidths();
        out.mwTransform(TopDown, false); // add coarse scale contributions
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
    bool absPrec;
    int maxScale;
};

