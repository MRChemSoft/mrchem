/** OperatorState is a simple helper class for operator application.
 * It keeps track of various state dependent variables and memory
 * regions. We cannot have some of this information directly in OperatorFunc
 * because of multi-threading issues.
*/

#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "MWNode.h"
#include "MathUtils.h"

#define GET_OP_IDX(FT,GT,ID) (2 * ((GT >> ID) & 1) + ((FT >> ID) & 1))

template<int D>
class OperatorState {
public:
    OperatorState(MWNode<D> &gn, double *scr1) : gNode(&gn) {
        this->kp1 = this->gNode->getKp1();
        this->kp1_d = this->gNode->getKp1_d();
        this->kp1_2 = MathUtils::ipow(this->kp1, 2);
        this->kp1_dm1 = MathUtils::ipow(this->kp1, D - 1);
        this->gData = this->gNode->getCoefs();
        this->maxDeltaL = -1;

        double *scr2 = scr1 + this->kp1_d;

        for (int i = 1; i < D; i++) {
            if (IS_ODD(i)) {
                this->aux[i] = scr2;
            } else {
                this->aux[i] = scr1;
            }
        }
    }
    void setFNode(MWNode<D> &fn) {
        this->fNode = &fn;
        this->fData = this->fNode->getCoefs();
        calcMaxDeltaL();
    }
    void setGComponent(int gt) {
        this->aux[D] = this->gData + gt * this->kp1_d;
        this->gt = gt;
    }
    void setFComponent(int ft) {
        this->aux[0] = this->fData + ft * this->kp1_d;
        this->ft = ft;
    }

    int getMaxDeltaL() const { return this->maxDeltaL; }
    int getOperIndex(int i) const { return GET_OP_IDX(this->ft, this->gt, i); }

    double **getAuxData() { return this->aux; }
    double **getOperData() { return this->oData; }

    friend class ConvolutionCalculator<D>;
    friend class DerivativeCalculator<D>;

private:
    int ft;
    int gt;
    int maxDeltaL;
    double fThreshold;
    double gThreshold;
    // Shorthands
    int kp1;
    int kp1_2;
    int kp1_d;
    int kp1_dm1;

    const OperatorTree *oTree;
    MWNode<D> *gNode;
    MWNode<D> *fNode;

    double *aux[D + 1];
    double *gData;
    double *fData;
    double *oData[D];

    void calcMaxDeltaL() {
        const int *gl = this->gNode->getNodeIndex().getTranslation();
        const int *fl = this->fNode->getNodeIndex().getTranslation();
        int max_dl = 0;
        for (int d = 0; d < D; d++) {
            int dl = abs(fl[d] - gl[d]);
            if (dl > max_dl) {
                max_dl = dl;
            }
        }
        this->maxDeltaL = max_dl;
    }
};

