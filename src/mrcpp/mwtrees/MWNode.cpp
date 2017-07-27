/**
 *  Simple n-dimensional node
 *
 *  Created on: Oct 12, 2009
 *      Author: jonas
 */

#include "MWNode.h"
#include "MWTree.h"
#include "ProjectedNode.h"
#include "MathUtils.h"
#include "QuadratureCache.h"
#include "GenNode.h"
#include "Timer.h"

using namespace std;
using namespace Eigen;

/** MWNode default constructor.
 *  Should be used only by SerialTree to obtain
 *  virtual table pointers for the derived classes. */
template<int D>
MWNode<D>::MWNode()
        : tree(0),
          parent(0),
          nodeIndex(),
          hilbertPath(),
          squareNorm(-1.0),
          n_coefs(0),
          coefs(0),
          status(0) {
    setIsLeafNode();
    setIsLooseNode();

    clearNorms();
    for (int i = 0; i < getTDim(); i++) {
        this->children[i] = 0;
    }

#ifdef HAVE_OPENMP
    omp_init_lock(&node_lock);
#endif
}

/** MWNode copy constructor.
 *  Creates loose nodes and copy coefs */
template<int D>
MWNode<D>::MWNode(const MWNode<D> &node)
        : tree(node.tree),
          parent(0),
          nodeIndex(node.nodeIndex),
          hilbertPath(node.hilbertPath),
          squareNorm(-1.0),
          n_coefs(0),
          coefs(0),
          status(0) {
    setIsLeafNode();
    setIsLooseNode();

    allocCoefs(this->getTDim(), this->getKp1_d());

    if (node.hasCoefs()) {
        setCoefBlock(0, node.getNCoefs(), node.getCoefs());
        if (this->getNCoefs() > node.getNCoefs()) {
            for(int i = node.getNCoefs(); i < this->getNCoefs(); i++) {
                this->coefs[i] = 0.0;
            }
        }
        this->setHasCoefs();
        this->calcNorms();
    } else {
        this->clearHasCoefs();
        this->clearNorms();
    }
    for (int i = 0; i < getTDim(); i++) {
        this->children[i] = 0;
    }

#ifdef HAVE_OPENMP
    omp_init_lock(&node_lock);
#endif
}

/** MWNode destructor.
  * Recursive deallocation of a node and all its decendants */
template<int D>
MWNode<D>::~MWNode() {
    if (this->isLooseNode()) this->freeCoefs();
#ifdef HAVE_OPENMP
    omp_destroy_lock(&node_lock);
#endif
}

/** Allocate the coefs vector. Only used by loose nodes. */
template<int D>
void MWNode<D>::allocCoefs(int n_blocks, int block_size) {
    if (this->n_coefs != 0) MSG_FATAL("n_coefs should be zero");
    if (this->isAllocated()) MSG_FATAL("Coefs already allocated");
    if (not this->isLooseNode()) MSG_FATAL("Only loose nodes here!");

    this->n_coefs = n_blocks * block_size;
    this->coefs = new double[this->n_coefs];

    this->clearHasCoefs();
    this->setIsAllocated();
}

/** Deallocation of coefficients. Only used by loose nodes. */
template<int D>
void MWNode<D>::freeCoefs() {
    if (not this->isLooseNode()) MSG_FATAL("Only loose nodes here!");

    if (this->coefs != 0) delete[] this->coefs;

    this->coefs = 0;
    this->n_coefs = 0;

    this->clearHasCoefs();
    this->clearIsAllocated();
}

template<int D>
void MWNode<D>::printCoefs() const {
    if (not this->isAllocated()) MSG_FATAL("Node is not allocated");
    println(0, "\nMW coefs");
    int kp1_d = this->getKp1_d();
    for (int i = 0; i < this->n_coefs; i++) {
        println(0, this->coefs[i]);
        if (i%kp1_d == 0) println(0, "\n");
    }
}

template<int D>
void MWNode<D>::getCoefs(Eigen::VectorXd &c) const {
    if (not this->isAllocated()) MSG_FATAL("Node is not allocated");
    if (not this->hasCoefs()) MSG_FATAL("Node has no coefs");
    if (this->n_coefs == 0) MSG_FATAL("ncoefs == 0");

    c = VectorXd::Map(this->coefs, this->n_coefs);
}

template<int D>
void MWNode<D>::zeroCoefs() {
    if (not this->isAllocated()) MSG_FATAL("Coefs not allocated");

    for (int i = 0; i < this->n_coefs; i++) {
        this->coefs[i] = 0.0;
    }
    this->zeroNorms();
    this->setHasCoefs();
}

template<int D>
void MWNode<D>::setCoefBlock(int block, int block_size, const double *c) {
    if (not this->isAllocated()) MSG_FATAL("Coefs not allocated");
    for (int i = 0; i < block_size; i++) {
        this->coefs[block*block_size + i] = c[i];
    }
}

template<int D>
void MWNode<D>::addCoefBlock(int block, int block_size, const double *c) {
    if (not this->isAllocated()) MSG_FATAL("Coefs not allocated");
    for (int i = 0; i < block_size; i++) {
        this->coefs[block*block_size + i] += c[i];
    }
}

template<int D>
void MWNode<D>::zeroCoefBlock(int block, int block_size) {
    if (not this->isAllocated()) MSG_FATAL("Coefs not allocated");
    for (int i = 0; i < block_size; i++) {
        this->coefs[block*block_size + i] = 0.0;
    }
}

template<int D>
void MWNode<D>::giveChildrenCoefs(bool overwrite) {
    assert(this->isBranchNode());
    if (not this->hasCoefs()) MSG_FATAL("No coefficients!");

    if (overwrite) {
        for (int i = 0; i < this->getTDim(); i++) {
            this->getMWChild(i).zeroCoefs();
        }
    }

    //coeff of child should be have been allocated already here
    int stride = this->getMWChild(0).getNCoefs();
    double* inp  = this->getCoefs();
    double* out = this->getMWChild(0).getCoefs();
    bool readOnlyScaling = false;
    if (this->isGenNode()) readOnlyScaling = true;

    this->tree->getSerialTree()->S_mwTransform(inp, out, readOnlyScaling, stride, overwrite);

    for (int i = 0; i < this->getTDim(); i++){
	this->getMWChild(i).setHasCoefs();
	this->getMWChild(i).calcNorms();//should need to compute only scaling norms
    }
}

/** Takes the scaling coefficients of the children and stores them consecutively
  * in the  given vector. */
template<int D>
void MWNode<D>::copyCoefsFromChildren() {
    int kp1_d = this->getKp1_d();
    int nChildren = this->getTDim();
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        MWNode<D> &child = getMWChild(cIdx);
        if (not child.hasCoefs()) MSG_FATAL("Child has no coefs");
        setCoefBlock(cIdx, kp1_d, child.getCoefs());
    }
}


/** Coefficient-Value transform
  *
  * This routine transforms the scaling coefficients of the node to the
  * function values in the corresponding quadrature roots (of its children).
  * Input parameter = forward: coef->value.
  * Input parameter = backward: value->coef.
  *
  * NOTE: this routine assumes a 0/1 (scaling on children 0 and 1)
  *       representation, in oppose to s/d (scaling and wavelet). */
template<int D>
void MWNode<D>::cvTransform(int operation) {
    int kp1 = this->getKp1();
    int kp1_dm1 = MathUtils::ipow(kp1, D - 1);
    int kp1_d = this->getKp1_d();
    int nCoefs = this->getTDim()*kp1_d;

    const ScalingBasis &sf = this->getMWTree().getMRA().getScalingBasis();
    const MatrixXd &S = sf.getCVMap(operation);

    double o_vec[nCoefs];
    double *out_vec = o_vec;
    double *in_vec = this->coefs;

    for (int i = 0; i < D; i++) {
        for (int t = 0; t < this->getTDim(); t++) {
            double *out = out_vec + t*kp1_d;
            double *in = in_vec + t*kp1_d;
            MathUtils::applyFilter(out, in, S, kp1, kp1_dm1, 0.0);
        }
        double *tmp = in_vec;
        in_vec = out_vec;
        out_vec = tmp;
    }
    int np1 = getScale() + 1; // we're working on scaling coefs on next scale
    double two_fac = pow(2.0, D*np1);
    if (operation == Backward) {
        two_fac = sqrt(1.0/two_fac);
    } else {
        two_fac = sqrt(two_fac);
    }
    if (IS_ODD(D)) {
        for (int i = 0; i < nCoefs; i++) {
            this->coefs[i] = two_fac*in_vec[i];
        }
    } else {
        for (int i = 0; i < nCoefs; i++) {
            this->coefs[i] *= two_fac;
        }
    }
}
/* Old interpolating version, somewhat faster
template<int D>
void MWNode<D>::cvTransform(int operation) {
    const ScalingBasis &sf = this->getMWTree().getMRA().getScalingBasis();
    if (sf.getScalingType() != Interpol) {
        NOT_IMPLEMENTED_ABORT;
    }

    int quadratureOrder = sf.getQuadratureOrder();
    getQuadratureCache(qc);

    double two_scale = pow(2.0, this->getScale() + 1);
    VectorXd modWeights = qc.getWeights(quadratureOrder);
    if (operation == Forward) {
        modWeights = modWeights.array().inverse();
        modWeights *= two_scale;
        modWeights = modWeights.array().sqrt();
    } else if (operation == Backward) {
        modWeights *= 1.0/two_scale;
        modWeights = modWeights.array().sqrt();
    } else {
        MSG_FATAL("Invalid operation");
    }

    int kp1 = this->getKp1();
    int kp1_d = this->getKp1_d();
    int kp1_p[D];
    for (int d = 0; d < D; d++) {
        kp1_p[d] = MathUtils::ipow(kp1, d);
    }

    for (int m = 0; m < this->getTDim(); m++) {
        for (int p = 0; p < D; p++) {
            int n = 0;
            for (int i = 0; i < kp1_p[D - p - 1]; i++) {
                for (int j = 0; j < kp1; j++) {
                    for (int k = 0; k < kp1_p[p]; k++) {
                        this->coefs[m * kp1_d + n] *= modWeights[j];
                        n++;
                    }
                }
            }
        }
    }
}
*/

/** Multiwavelet transform: fast version
  *
  * Application of the filters on one node to pass from a 0/1 (scaling
  * on children 0 and 1) representation to an s/d (scaling and
  * wavelet) representation. Bit manipulation is used in order to
  * determine the correct filters and whether to apply them or just
  * pass to the next couple of indexes. The starting coefficients are
  * preserved until the application is terminated, then they are
  * overwritten. With minor modifications this code can also be used
  * for the inverse mw transform (just use the transpose filters) or
  * for the application of an operator (using A, B, C and T parts of an
  * operator instead of G1, G0, H1, H0). This is the version where the
  * three directions are operated one after the other. Although this
  * is formally faster than the other algorithm, the separation of the
  * three dimensions prevent the possibility to use the norm of the
  * operator in order to discard a priori negligible contributions.

  * Luca Frediani, August 2006
  * C++ version: Jonas Juselius, September 2009 */
template<int D>
void MWNode<D>::mwTransform(int operation) {
    int kp1 = this->getKp1();
    int kp1_dm1 = MathUtils::ipow(kp1, D - 1);
    int kp1_d = this->getKp1_d();
    int nCoefs = this->getTDim()*kp1_d;
    const MWFilter &filter = getMWTree().getMRA().getFilter();
    double overwrite = 0.0;

    double o_vec[nCoefs];
    double *out_vec = o_vec;
    double *in_vec = this->coefs;

    for (int i = 0; i < D; i++) {
        int mask = 1 << i;
        for (int gt = 0; gt < this->getTDim(); gt++) {
            double *out = out_vec + gt * kp1_d;
            for (int ft = 0; ft < this->getTDim(); ft++) {
                /* Operate in direction i only if the bits along other
                 * directions are identical. The bit of the direction we
                 * operate on determines the appropriate filter/operator */
                if ((gt | mask) == (ft | mask)) {
                    double *in = in_vec + ft * kp1_d;
                    int fIdx = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                    const MatrixXd &oper = filter.getSubFilter(fIdx, operation);
                    MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
                    overwrite = 1.0;
                }
            }
            overwrite = 0.0;
        }
        double *tmp = in_vec;
        in_vec = out_vec;
        out_vec = tmp;
    }
    if (IS_ODD(D)) {
        for (int i = 0; i < nCoefs; i++) {
            this->coefs[i] = in_vec[i];
        }
    }
}

/** Set all norms to Undefined. */
template<int D>
void MWNode<D>::clearNorms() {
    this->squareNorm = -1.0;
    for (int i = 0; i < this->getTDim(); i++) {
        this->componentNorms[i] = -1.0;
    }
}

/** Set all norms to zero. */
template<int D>
void MWNode<D>::zeroNorms() {
    this->squareNorm = 0.0;
    for (int i = 0; i < this->getTDim(); i++) {
        this->componentNorms[i] = 0.0;
    }
}

/** Calculate and store square norm and component norms, if allocated. */
template<int D>
void MWNode<D>::calcNorms() {
    this->squareNorm = 0.0;
    for (int i = 0; i < this->getTDim(); i++) {
        double norm_i = calcComponentNorm(i);
        this->componentNorms[i] = norm_i;
        this->squareNorm += norm_i*norm_i;
    }
}

/** Calculate and return the squared scaling norm. */
template<int D>
double MWNode<D>::getScalingNorm() const {
    double sNorm = this->getComponentNorm(0);
    if (sNorm >= 0.0) {
        return sNorm*sNorm;
    } else {
        return -1.0;
    }
}

/** Calculate and return the squared wavelet norm. */
template<int D>
double MWNode<D>::getWaveletNorm() const {
    double wNorm = 0.0;
    for (int i = 1; i < this->getTDim(); i++) {
        double norm_i = this->getComponentNorm(i);
        if (norm_i >= 0.0) {
            wNorm += norm_i*norm_i;
        } else {
            wNorm = -1.0;
        }
    }
    return wNorm;
}

/** Calculate the norm of one component (NOT the squared norm!). */
template<int D>
double MWNode<D>::calcComponentNorm(int i) const {
  assert(i==0 or (not this->isGenNode()));//GenNodes have no wavelets coefficients
    assert(this->isAllocated());
    assert(this->hasCoefs());

    const double *c = this->getCoefs();
    int size = this->getKp1_d();
    int start = i*size;

    double sq_norm = 0.0;
#ifdef HAVE_BLAS
    sq_norm = cblas_ddot(size, &c[start], 1, &c[start], 1);
#else
    for (int i = start; i < start+size; i++) {
        sq_norm += c[i]*c[i];
    }
#endif
    return sqrt(sq_norm);
}

/** Update the coefficients of the node by a mw transform of the scaling
  * coefficients of the children. Option to overwrite or add up existing
  * coefficients. */
template<int D>
void MWNode<D>::reCompress() {
    if (this->isBranchNode()) {
        if (not this->isAllocated()) MSG_FATAL("Coefs not allocated");
        copyCoefsFromChildren();
        mwTransform(Compression);
        this->setHasCoefs();
        this->calcNorms();
    }
}

/** Recurse down until an EndNode is found, and then crop children with
  * too high precision. */
template<int D>
bool MWNode<D>::crop(double prec, double splitFac, bool absPrec) {
    if (this->isEndNode()) {
        return true;
    } else {
        for (int i = 0; i < this->getTDim(); i++) {
            MWNode<D> &child = *this->children[i];
            if (child.crop(prec, splitFac, absPrec)) {
                if (this->splitCheck(prec, splitFac, absPrec) == false) {
                    this->deleteChildren();
                    return true;
                }
            }
        }
    }
    return false;
}

template<int D>
bool MWNode<D>::splitCheck(double prec, double splitFac, bool absPrec) const {
    if (prec < 0.0) {
        return false;
    }
    double scale_fac = getScaleFactor(splitFac, absPrec);
    double w_thrs = max(2.0*MachinePrec, prec*scale_fac);
    double w_norm = sqrt(getWaveletNorm());
    if (w_norm > w_thrs) {
        return true;
    } else {
        return false;
    }
}

/** Calculate the threshold for the wavelet norm.
  *
  * Calculates the threshold that has to be met in the wavelet norm in order to
  * guarantee the precision in the function representation. Depends on the
  * square norm of the function and the requested relative accuracy. */
template<int D>
double MWNode<D>::getScaleFactor(double splitFac, bool absPrec) const {
    double t_norm = 1.0;
    double sq_norm = this->tree->getSquareNorm();
    if (sq_norm > 0.0 and not absPrec) {
        t_norm = sqrt(sq_norm);
    }
    double scale_fac = 1.0;
    if (splitFac > MachineZero) {
        double expo = 0.5 * splitFac * (getScale() + 1);
        scale_fac = pow(2.0, -expo);
    }
    return t_norm * scale_fac;
}

template<int D>
void MWNode<D>::createChildren() {
    if (this->isBranchNode()) MSG_FATAL("Node already has children");
    this->getMWTree().getSerialTree()->allocChildren(*this);
    this->setIsBranchNode();
}

/** Recursive deallocation of children and all their descendants.
  * Leaves node as LeafNode and children[] as null pointer. */
template<int D>
void MWNode<D>::deleteChildren() {
    if (this->isLeafNode()) return;
    for (int cIdx = 0; cIdx < getTDim(); cIdx++) {
        if (this->children[cIdx] != 0) {
	    MWNode<D> &child = getMWChild(cIdx);
	    child.deleteChildren();
            child.dealloc();
	}
	this->children[cIdx] = 0;
    }
    this->setIsLeafNode();
}

template<int D>
void MWNode<D>::deleteGenerated() {
    if (this->isBranchNode()) {      
        if (this->isEndNode()) {
            this->deleteChildren();
        } else {
            for (int cIdx = 0; cIdx < getTDim(); cIdx++) {
                this->getMWChild(cIdx).deleteGenerated();
            }
        }
    }
}

template<int D>
void MWNode<D>::getCenter(double *r) const {
    NOT_IMPLEMENTED_ABORT;
    //    assert(r != 0);
    //    double sFac = pow(2.0, -getScale());
    //    for (int d = 0; d < D; d++) {
    //        double l = (double) getTranslation()[d];
    //        r[d] = sFac*(l + 0.5);
    //    }
}

template<int D>
void MWNode<D>::getBounds(double *lb, double *ub) const {
    int n = getScale();
    double p = pow(2.0, -n);
    const int *l = getTranslation();
    for (int i = 0; i < D; i++) {
        lb[i] = p * l[i];
        ub[i] = p * (l[i] + 1);
    }
}

/** Routine to find the path along the tree.
  *
  * Given the translation indices at the final scale, computes the child m
  * to be followed at the current scale in oder to get to the requested
  * node at the final scale. The result is the index of the child needed.
  * The index is obtained by bit manipulation of of the translation indices. */
template<int D>
int MWNode<D>::getChildIndex(const NodeIndex<D> &nIdx) const {
    assert(isAncestor(nIdx));
    int cIdx = 0;
    int diffScale = nIdx.getScale() - getScale() - 1;
    assert(diffScale >= 0);
    for (int d = 0; d < D; d++) {
        int bit = (nIdx.getTranslation()[d] >> (diffScale)) & 1;
        cIdx = cIdx + (bit << d);
    }
    assert(cIdx >= 0);
    assert(cIdx < getTDim());
    return cIdx;
}

/** Routine to find the path along the tree.
  *
  * Given a point in space, determines which child should be followed
  * to get to the corresponding terminal node. */
template<int D>
int MWNode<D>::getChildIndex(const double *r) const {
    assert(hasCoord(r));
    int cIdx = 0;
    double sFac = pow(2.0, -getScale());
    const int *l = getTranslation();
    for (int d = 0; d < D; d++) {
        if (r[d] > sFac*(l[d] + 0.5)) {
            cIdx = cIdx + (1 << d);
        }
    }
    assert(cIdx >= 0);
    assert(cIdx < getTDim());
    return cIdx;
}

template<int D>
void MWNode<D>::getPrimitiveQuadPts(MatrixXd &pts) const {
    int kp1 = this->getKp1();
    pts = MatrixXd::Zero(kp1,D);

    getQuadratureCache(qc);
    const VectorXd &roots = qc.getRoots(kp1);

    double sFac = pow(2.0, -this->getScale());
    const int *l = this->getTranslation();
    for (int d = 0; d < D; d++) {
        pts.col(d) = sFac*(roots.array() + double(l[d]));
    }
}

template<int D>
void MWNode<D>::getPrimitiveChildPts(MatrixXd &pts) const {
    int kp1 = this->getKp1();
    pts = MatrixXd::Zero(D,2*kp1);

    getQuadratureCache(qc);
    const VectorXd &roots = qc.getRoots(kp1);

    double sFac = pow(2.0, -(this->getScale() + 1));
    const int *l = this->getTranslation();
    for (int d = 0; d < D; d++) {
        pts.row(d).segment(0, kp1) = sFac*(roots.array() + 2.0*double(l[d]));
        pts.row(d).segment(kp1, kp1) = sFac*(roots.array() + 2.0*double(l[d]) + 1);
    }
}

template<int D>
void MWNode<D>::getExpandedQuadPts(Eigen::MatrixXd &pts) const {
    MatrixXd prim_pts;
    getPrimitiveQuadPts(prim_pts);

    int kp1 = this->getKp1();
    int kp1_d = this->getKp1_d();
    pts = MatrixXd::Zero(D, kp1_d);

    if (D == 1) pts = prim_pts;
    if (D == 2) MathUtils::tensorExpandCoords_2D(kp1, prim_pts, pts);
    if (D == 3) MathUtils::tensorExpandCoords_3D(kp1, prim_pts, pts);
    if (D >= 4) NOT_IMPLEMENTED_ABORT;
}

template<int D>
void MWNode<D>::getExpandedChildPts(MatrixXd &pts) const {
    MatrixXd prim_pts;
    getPrimitiveChildPts(prim_pts);

    int tDim = this->getTDim();
    int kp1 = this->getKp1();
    int kp1_d = this->getKp1_d();
    pts = MatrixXd::Zero(D, tDim*kp1_d);
    MatrixXd prim_t = MatrixXd::Zero(D, kp1);
    MatrixXd exp_t = MatrixXd::Zero(D, kp1_d);

    for (int t = 0; t < tDim; t++) {
        for (int d = 0; d < D; d++) {
            int idx = (t>>d)&1;
            prim_t.row(d) = prim_pts.block(d, idx*kp1, 1, kp1);
        }
        if (D == 1) exp_t = prim_t;
        if (D == 2) MathUtils::tensorExpandCoords_2D(kp1, prim_t, exp_t);
        if (D == 3) MathUtils::tensorExpandCoords_3D(kp1, prim_t, exp_t);
        if (D >= 4) NOT_IMPLEMENTED_ABORT;
        pts.block(0, t*kp1_d, D, kp1_d) = exp_t;
    }
}

/** Const version of node retriever that NEVER generates.
  *
  * Recursive routine to find and return the node with a given NodeIndex.
  * This routine returns the appropriate ProjectedNode, or a NULL pointer if
  * the node does not exist, or if it is a GenNode. Recursion starts at at this
  * node and ASSUMES the requested node is in fact decending from this node. */
template<int D>
const MWNode<D> *MWNode<D>::retrieveNodeNoGen(const NodeIndex<D> &idx) const {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(this->isAncestor(idx));
    if (this->isEndNode()) { // don't return GenNodes
        return 0;
    }
    int cIdx = getChildIndex(idx);
    assert(this->children[cIdx] != 0);
    return this->children[cIdx]->retrieveNodeNoGen(idx);
}

/** Node retriever that NEVER generates.
  *
  * Recursive routine to find and return the node with a given NodeIndex.
  * This routine returns the appropriate ProjectedNode, or a NULL pointer if
  * the node does not exist, or if it is a GenNode. Recursion starts at at this
  * node and ASSUMES the requested node is in fact decending from this node. */
template<int D>
MWNode<D> *MWNode<D>::retrieveNodeNoGen(const NodeIndex<D> &idx) {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(this->isAncestor(idx));
    if (this->isEndNode()) { // don't return GenNodes
        return 0;
    }
    int cIdx = getChildIndex(idx);
    assert(this->children[cIdx] != 0);
    return this->children[cIdx]->retrieveNodeNoGen(idx);
}

template<int D>
const MWNode<D> *MWNode<D>::retrieveNodeOrEndNode(const double *r, int depth) const {
    if (getDepth() == depth or this->isEndNode()) {
        return this;
    }
    int cIdx = getChildIndex(r);
    assert(this->children[cIdx] != 0);
    return this->children[cIdx]->retrieveNodeOrEndNode(r, depth);
}

/** Node retriever that return requested ProjectedNode or EndNode.
  *
  * Recursive routine to find and return the node with a given NodeIndex.
  * This routine returns the appropriate ProjectedNode, or the EndNode on the
  * path to the requested node, and will never create or return GenNodes.
  * Recursion starts at at this node and ASSUMES the requested node is in fact
  * decending from this node. */
template<int D>
MWNode<D> *MWNode<D>::retrieveNodeOrEndNode(const double *r, int depth) {
    if (getDepth() == depth or this->isEndNode()) {
        return this;
    }
    int cIdx = getChildIndex(r);
    assert(this->children[cIdx] != 0);
    return this->children[cIdx]->retrieveNodeOrEndNode(r, depth);
}

template<int D>
const MWNode<D> *MWNode<D>::retrieveNodeOrEndNode(const NodeIndex<D> &idx) const {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(isAncestor(idx));
    // We should in principle lock before read, but it makes things slower,
    // and the EndNode status does not change (normally ;)
    if (isEndNode()) {
        return this;
    }
    int cIdx = getChildIndex(idx);
    assert(children[cIdx] != 0);
    return this->children[cIdx]->retrieveNodeOrEndNode(idx);
}

template<int D>
MWNode<D> *MWNode<D>::retrieveNodeOrEndNode(const NodeIndex<D> &idx) {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(isAncestor(idx));
    // We should in principle lock before read, but it makes things slower,
    // and the EndNode status does not change (normally ;)
    if (isEndNode()) {
        return this;
    }
    int cIdx = getChildIndex(idx);
    assert(children[cIdx] != 0);
    return this->children[cIdx]->retrieveNodeOrEndNode(idx);
}

/** Node retriever that ALWAYS returns the requested node.
  *
  * Recursive routine to find and return the node with a given NodeIndex.
  * This routine always returns the appropriate node, and will generate nodes
  * that does not exist. Recursion starts at this node and ASSUMES the
  * requested node is in fact decending from this node. */
template<int D>
MWNode<D> *MWNode<D>::retrieveNode(const double *r, int depth) {
    if (depth < 0) MSG_FATAL("Invalid argument");

    if (getDepth() == depth) {
        return this;
    }
    assert(hasCoord(r));
    // If we have reached an endNode, lock if necessary, and start generating
    // NB! retrieveNode() for GenNodes behave a bit differently.
    SET_NODE_LOCK();
    if (this->isLeafNode()) {
        genChildren();
        giveChildrenCoefs();
   }
    UNSET_NODE_LOCK();
    int cIdx = getChildIndex(r);
    assert(this->children[cIdx] != 0);
    return this->children[cIdx]->retrieveNode(r, depth);
}

/** Node retriever that ALWAYS returns the requested node, possibly without coefs.
  *
  * Recursive routine to find and return the node with a given NodeIndex. This
  * routine always returns the appropriate node, and will generate nodes that
  * does not exist. Recursion starts at this node and ASSUMES the requested
  * node is in fact decending from this node. */
template<int D>
MWNode<D> *MWNode<D>::retrieveNode(const NodeIndex<D> &idx) {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(isAncestor(idx));
    SET_NODE_LOCK();
    if (isLeafNode()) {
      genChildren();
      giveChildrenCoefs();
    }
    UNSET_NODE_LOCK();
    int cIdx = getChildIndex(idx);

    assert(this->children[cIdx] != 0);
    return this->children[cIdx]->retrieveNode(idx);
}

/** Test if a given coordinate is within the boundaries of the node. */
template<int D>
bool MWNode<D>::hasCoord(const double *r) const {
    double sFac = pow(2.0, -getScale());
    const int *l = getTranslation();
    //    println(1, "[" << r[0] << "," << r[1] << "," << r[2] << "]");
    //    println(1, "[" << l[0] << "," << l[1] << "," << l[2] << "]");
    //    println(1, *this);
    for (int d = 0; d < D; d++) {
        if (r[d] < sFac*l[d] or r[d] > sFac*(l[d] + 1)) {
            //            println(1, "false");
            return false;
        }
    }
    //    println(1, "true");
    return true;
}

/** Testing if nodes are compatible wrt NodeIndex and Tree (order, rootScale,
  * relPrec, etc). */
template<int D>
bool MWNode<D>::isCompatible(const MWNode<D> &node) {
    NOT_IMPLEMENTED_ABORT;
    //    if (nodeIndex != node.nodeIndex) {
    //        println(0, "nodeIndex mismatch" << std::endl);
    //        return false;
    //    }
    //    if (not this->tree->checkCompatible(*node.tree)) {
    //        println(0, "tree type mismatch" << std::endl);
    //        return false;
    //    }
    //    return true;
}

/** Test if the node is decending from a given NodeIndex, that is, if they have
  * overlapping support. */
template<int D>
bool MWNode<D>::isAncestor(const NodeIndex<D> &idx) const {
    int relScale = idx.getScale() - getScale();
    if (relScale < 0) {
        return false;
    }
    const int *l = getTranslation();
    for (int d = 0; d < D; d++) {
        int reqTransl = idx.getTranslation()[d] >> relScale;
        if (l[d] != reqTransl) {
            return false;
        }
    }
    return true;
}

template<int D>
bool MWNode<D>::isDecendant(const NodeIndex<D> &idx) const {
    NOT_IMPLEMENTED_ABORT;
}

template class MWNode<1>;
template class MWNode<2>;
template class MWNode<3>;
