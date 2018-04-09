#include "DerivativeCalculator.h"
#include "DerivativeOperator.h"
#include "OperatorState.h"
#include "OperatorNode.h"
#include "FunctionTree.h"
#include "BandWidth.h"
#include "Timer.h"
#include "eigen_disable_warnings.h"
#include "blas.h"

using namespace std;
using namespace Eigen;

template<int D>
DerivativeCalculator<D>::DerivativeCalculator(int dir,
                                              DerivativeOperator<D> &o,
                                              FunctionTree<D> &f)
        : applyDir(dir),
          fTree(&f),
          oper(&o) {
    if (dir < 0 or dir >= D) MSG_FATAL("Invalid apply dir");
    initTimers();
}

template<int D>
DerivativeCalculator<D>::~DerivativeCalculator() {
    this->operStat.flushNodeCounters();
    println(10, this->operStat);
}

template<int D>
void DerivativeCalculator<D>::initTimers() {
    int nThreads = omp_get_max_threads();
    for (int i = 0; i < nThreads; i++) {
        this->band_t.push_back(Timer(false));
        this->calc_t.push_back(Timer(false));
        this->norm_t.push_back(Timer(false));
    }
}

template<int D>
void DerivativeCalculator<D>::clearTimers() {
    this->band_t.clear();
    this->calc_t.clear();
    this->norm_t.clear();
}

template<int D>
void DerivativeCalculator<D>::printTimers() const {
    int oldprec = TelePrompter::setPrecision(1);
    int nThreads = omp_get_max_threads();
    printout(20, endl);
    printout(20, endl << "thread ");
    for (int i = 0; i < nThreads; i++) {
        printout(20, setw(9) << i);
    }
    printout(20, endl << "band     ");
    for (int i = 0; i < nThreads; i++) {
        printout(20, this->band_t[i].getWallTime() << "  ");
    }
    printout(20, endl << "calc     ");
    for (int i = 0; i < nThreads; i++) {
        printout(20, this->calc_t[i].getWallTime() << "  ");
    }
    printout(20, endl << "norm     ");
    for (int i = 0; i < nThreads; i++) {
        printout(20, this->norm_t[i].getWallTime() << "  ");
    }
    printout(20, endl);
    printout(20, endl);
    TelePrompter::setPrecision(oldprec);
}

template<int D>
void DerivativeCalculator<D>::calcNode(MWNode<D> &gNode) {
    gNode.zeroCoefs();

    int nComp = (1 << D);
    int depth = gNode.getDepth();
    double tmpCoefs[gNode.getNCoefs()];
    OperatorState<D> os(gNode, tmpCoefs);
    this->operStat.incrementGNodeCounters(gNode);

    // Get all nodes in f within the bandwith of O in g
    this->band_t[omp_get_thread_num()].resume();
    MWNodeVector fBand = makeOperBand(gNode);
    this->band_t[omp_get_thread_num()].stop();

    assert(this->oper->size() == 1);
    const OperatorTree &oTree = this->oper->getComponent(0);
    os.oTree = &oTree;

    this->calc_t[omp_get_thread_num()].resume();
    for (int n = 0; n < fBand.size(); n++) {
        MWNode<D> &fNode = *fBand[n];
        os.setFNode(fNode);
        for (int ft = 0; ft < nComp; ft++) {
            double fNorm = fNode.getComponentNorm(ft);
            if (fNorm < MachineZero) {
                continue;
            }
            os.setFComponent(ft);
            for (int gt = 0; gt < nComp; gt++) {
                os.setGComponent(gt);
                applyOperator(os);
            }
        }
    }
    this->calc_t[omp_get_thread_num()].stop();

    this->norm_t[omp_get_thread_num()].resume();
    gNode.calcNorms();
    this->norm_t[omp_get_thread_num()].stop();
}

/** Return a vector of nodes in F affected by O, given a node in G */
template<int D>
MWNodeVector DerivativeCalculator<D>::makeOperBand(const MWNode<D> &gNode) {
    assert(this->applyDir >= 0);
    assert(this->applyDir < D);

    MWNodeVector band;
    const NodeIndex<D> &idx_0 = gNode.getNodeIndex();

    // Assumes given width only in applyDir, otherwise width = 0
    int width = this->oper->getMaxBandWidth();
    for (int w = -width; w <= width; w++) { 
        NodeIndex<D> idx_w(idx_0);
        idx_w.getTranslation()[this->applyDir] += w;

        // returns -1 if out of bounds
        int rIdx_w = this->fTree->getRootIndex(idx_w);
        if (rIdx_w >= 0) band.push_back(&this->fTree->getNode(idx_w));
    }
    return band;
}

/** Apply a single operator component (term) to a single f-node. Whether the
operator actualy is applied is determined by a screening threshold. */
template<int D>
void DerivativeCalculator<D>::applyOperator(OperatorState<D> &os) {
    const OperatorTree &oTree = *os.oTree;
    MWNode<D> &gNode = *os.gNode;
    MWNode<D> &fNode = *os.fNode;

    const int *gTransl = gNode.getTranslation();
    const int *fTransl = fNode.getTranslation();
    int depth = gNode.getDepth();

    double oNorm = 1.0;
    double **oData = os.getOperData();

    for (int d = 0; d < D; d++) {
        int oTransl = fTransl[d] - gTransl[d];

//  The following will check the actual band width in each direction.
//  Not needed if the thresholding at the end of this routine is active.
        int a = (os.gt & (1 << d)) >> d;
        int b = (os.ft & (1 << d)) >> d;
        int idx = (a << 1) + b;
        int w = oTree.getBandWidth().getWidth(depth, idx);
        if (abs(oTransl) > w) {
            return;
        }

        const OperatorNode &oNode = oTree.getNode(depth, oTransl);
        int oIdx = os.getOperIndex(d);
        double ocn = oNode.getComponentNorm(oIdx);
        oNorm *= ocn;
        if (this->applyDir == d) {
            oData[d] = const_cast<double *>(oNode.getCoefs()) + oIdx*os.kp1_2;
        } else {
            if (oTransl == 0 and (oIdx == 0 or oIdx == 3)) {
                // This will activate the identity operator in direction i
                oData[d] = 0;
            } else {
                // This means that we are in a zero part of the identity operator
                return;
            }
        }
    }
    this->operStat.incrementFNodeCounters(fNode, os.ft, os.gt);
    tensorApplyOperComp(os);
}

/** Perorm the required linear algebra operations in order to apply an
operator component to a f-node in a n-dimensional tesor space. */
template<int D>
void DerivativeCalculator<D>::tensorApplyOperComp(OperatorState<D> &os) {
    double **aux = os.getAuxData();
    double **oData = os.getOperData();
#ifdef HAVE_BLAS
    double mult = 0.0;
    for (int i = 0; i < D; i++) {
        if (oData[i] != 0) {
            if (i == D - 1) { // Last dir: Add up into g
                mult = 1.0;
            }
            const double *f = aux[i];
            double *g = const_cast<double *>(aux[i + 1]);
            dgemm('t', 'n', os.kp1_dm1, os.kp1, os.kp1, 1.0, f,
                  os.kp1, oData[i], os.kp1, mult, g, os.kp1_dm1);
        } else {
            // Identity operator in direction i
            Map<MatrixXd> f(aux[i], os.kp1, os.kp1_dm1);
            Map<MatrixXd> g(aux[i + 1], os.kp1_dm1, os.kp1);
            if (oData[i] == 0) {
                if (i == D - 1) { // Last dir: Add up into g
                    g += f.transpose();
                } else {
                    g = f.transpose();
                }
            }
        }
    }
#else
    for (int i = 0; i < D; i++) {
        Map<MatrixXd> f(aux[i], os.kp1, os.kp1_dm1);
        Map<MatrixXd> g(aux[i + 1], os.kp1_dm1, os.kp1);
        if (oData[i] != 0) {
            Eigen::Map<MatrixXd> op(oData[i], os.kp1, os.kp1);
            if (i == D - 1) { // Last dir: Add up into g
                g += f.transpose() * op;
            } else {
                g = f.transpose() * op;
            }
        } else {
            // Identity operator in direction i
            if (i == D - 1) { // Last dir: Add up into g
                g += f.transpose();
            } else {
                g = f.transpose();
            }
        }
    }
#endif
}

template<int D>
MWNodeVector* DerivativeCalculator<D>::getInitialWorkVector(MWTree<D> &tree) const {
    return tree.copyEndNodeTable();
}

template class DerivativeCalculator<1>;
template class DerivativeCalculator<2>;
template class DerivativeCalculator<3>;

