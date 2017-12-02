#pragma once

#include "TreeCalculator.h"
#include "OperatorStatistics.h"
#include "mrcpp_declarations.h"

template<int D>
class ConvolutionCalculator : public TreeCalculator<D> {
public:
    ConvolutionCalculator(double p,
                          ConvolutionOperator<D> &o,
                          FunctionTree<D> &f,
                          int depth = MaxDepth);
    virtual ~ConvolutionCalculator();

    virtual MWNodeVector* getInitialWorkVector(MWTree<D> &tree) const;

protected:
    int maxDepth;
    double prec;
    ConvolutionOperator<D> *oper;
    FunctionTree<D> *fTree;
    std::vector<Timer> band_t;
    std::vector<Timer> calc_t;
    std::vector<Timer> norm_t;

    OperatorStatistics<D> operStat;
    std::vector<Eigen::MatrixXi *> bandSizes;

    static const int nComp = (1 << D);
    static const int nComp2 = (1 << D) * (1 << D);

    MWNodeVector* makeOperBand(const MWNode<D> &gNode);
    void fillOperBand(MWNodeVector *band, NodeIndex<D> &idx, const int *nbox, int dim);

    void initTimers();
    void clearTimers();
    void printTimers() const;

    void initBandSizes();
    int getBandSizeFactor(int i, int depth,const OperatorState<D> &os) const;
    void calcBandSizeFactor(Eigen::MatrixXi &bs, int depth, const BandWidth &bw);

    virtual void calcNode(MWNode<D> &node);
    virtual void postProcess() {
        printTimers();
        clearTimers();
        initTimers();
    }

    void applyOperComp(OperatorState<D> &os);
    void applyOperator(OperatorState<D> &os);
    void tensorApplyOperComp(OperatorState<D> &os);
};


