#pragma once

#include "TreeBuilder.h"
#include "MultiplicationCalculator.h"
#include "WaveletAdaptor.h"
#include "Timer.h"

template<int D>
class MWMultiplier {
public:
    MWMultiplier(double pr = -1.0, int ms = MaxScale)
        : prec(pr), maxScale(ms) { }
    virtual ~MWMultiplier() { }

    double getPrecision() const { return this->prec; }
    int getMaxScale() const { return this->maxScale; }

    void setPrecision(double pr) { this->prec = pr; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    void operator()(FunctionTree<3> &out, double c,
                    FunctionTree<D> &tree_a,
                    FunctionTree<D> &tree_b,
                    int maxIter = -1) const {
        FunctionTreeVector<D> tree_vec;
        tree_vec.push_back(c, &tree_a);
        tree_vec.push_back(1.0, &tree_b);
        return (*this)(out, tree_vec, maxIter);
    }
    void operator()(FunctionTree<D> &out,
                    FunctionTreeVector<D> &inp,
                    int maxIter = -1) const {
        TreeBuilder<D> builder;
        WaveletAdaptor<D> adaptor(this->prec, this->maxScale);
        MultiplicationCalculator<D> calculator(inp);

        builder.build(out, calculator, adaptor, maxIter);

        Timer trans_t;
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        trans_t.stop();

        Timer clean_t;
        for (int i = 0; i < inp.size(); i++) {
            FunctionTree<D> &tree = inp.getFunc(i);
            tree.deleteGenerated();
        }
        clean_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, "Time cleaning       " << clean_t);
        println(10, std::endl);
    }
protected:
    double prec;
    int maxScale;
};

