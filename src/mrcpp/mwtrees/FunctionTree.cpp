/**
 *  \date 2016
 *          CTCC, University of Tromsø
 */

#include "FunctionTree.h"
#include "SerialFunctionTree.h"
#include "FunctionNode.h"
#include "ProjectedNode.h"
#include "HilbertIterator.h"

using namespace std;
using namespace Eigen;

/** FunctionTree constructor for Serial Tree.
  * */
template<int D>
FunctionTree<D>::FunctionTree(const MultiResolutionAnalysis<D> &mra)
        : MWTree<D> (mra) {
    this->serialTree_p = new SerialFunctionTree<D>(this);
    this->serialTree_p->allocRoots(*this);
    this->resetEndNodeTable();
}

/** FunctionTree constructor for Serial Tree using shared memory.
  * */
template<int D>
FunctionTree<D>::FunctionTree(const MultiResolutionAnalysis<D> &mra, SharedMemory* &shMem)
        : MWTree<D> (mra) {
    this->serialTree_p = new SerialFunctionTree<D>(this);
    this->serialTree_p->isShared = true;
    this->serialTree_p->shMem = shMem;    
    this->serialTree_p->allocRoots(*this);
    this->resetEndNodeTable();
}

//** FunctionTree destructor. */
template<int D>
FunctionTree<D>::~FunctionTree() {
    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<D> &root = this->getRootMWNode(i);
        root.deleteChildren();
        root.dealloc();
        this->rootBox.clearNode(i);
    }
    delete this->serialTree_p;
}

/** Leaves the tree inn the same state as after construction*/
template<int D>
void FunctionTree<D>::clear() {
    NOT_IMPLEMENTED_ABORT;
}

/** Write the tree structure to disk, for later use.
  * Argument file name will get a ".tree" file extension, and in MPI an
  * additional "-[rank]". */
template<int D>
bool FunctionTree<D>::saveTree(const string &file) {
    NOT_IMPLEMENTED_ABORT;
}

/** Read a previously stored tree structure from disk.
  * Argument file name will get a ".tree" file extension, and in MPI an
  * additional "-[rank]". */
template<int D>
bool FunctionTree<D>::loadTree(const string &file) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
double FunctionTree<D>::integrate() const {
    double result = 0.0;
    for (int i = 0; i < this->rootBox.size(); i++) {
        const FunctionNode<D> &fNode = getRootFuncNode(i);
        result += fNode.integrate();
    }
    return result;
}

template<int D>
double FunctionTree<D>::dot(const FunctionTree<D> &ket) {
    const FunctionTree<D> &bra = *this;
    if (bra.getMRA() != ket.getMRA()){
      cout<<ket.getMRA().getScalingBasis()<<" "<<bra.getMRA().getScalingBasis()<<endl;
      cout<<ket.getMRA().getWorldBox()<<" "<<bra.getMRA().getWorldBox()<<endl;
      cout<<ket.getMRA().getMaxDepth()<<" "<<bra.getMRA().getMaxDepth()<<endl;
      MSG_FATAL("Trees not compatible");
    }
    MWNodeVector nodeTable;
    HilbertIterator<D> it(this);
    it.setReturnGenNodes(false);
    while(it.next()) {
        MWNode<D> &node = it.getNode();
        nodeTable.push_back(&node);
    }
    int nNodes = nodeTable.size();
    double result = 0.0;
    double locResult = 0.0;
//OMP is disabled in order to get EXACT results (to the very last digit), the
//order of summation makes the result different beyond the 14th digit or so.
//OMP does improve the performace, but its not worth it for the time being.
//#pragma omp parallel firstprivate(n_nodes, locResult)
//		shared(nodeTable,rhs,result)
//    {
//#pragma omp for schedule(guided)
    for (int n = 0; n < nNodes; n++) {
        const FunctionNode<D> &braNode = static_cast<const FunctionNode<D> &>(*nodeTable[n]);
        const MWNode<D> *mwNode = ket.findNode(braNode.getNodeIndex());
        if (mwNode == 0) continue;

        const FunctionNode<D> &ketNode = static_cast<const FunctionNode<D> &>(*mwNode);
        if (braNode.isRootNode()) {
            locResult += braNode.dotScaling(ketNode);
        }
        locResult += braNode.dotWavelet(ketNode);
    }
//#pragma omp critical
    result += locResult;
//    }
    return result;
}

template<int D>
double FunctionTree<D>::evalf(const double *r) {
    MWNode<D> &mr_node = this->getNodeOrEndNode(r);
    FunctionNode<D> &f_node = static_cast<FunctionNode<D> &>(mr_node);
    double result = f_node.evalf(r);
    this->deleteGenerated();
    return result;
}

template<int D>
void FunctionTree<D>::square() {
    NOT_IMPLEMENTED_ABORT
    // Doesn't work in MPI
//    this->purgeGenNodes();
//    for (int i = 0; i < this->endNodeTable.size(); i++) {
//        MWNode<D> &node = *this->endNodeTable[i];
//        node.mwTransform(Reconstruction);
//        node.cvTransform(MWNode<D>::Forward);
//        node.getCoefs() = node.getCoefs().array().square();
//        node.cvTransform(MWNode<D>::Backward);
//        node.mwTransform(Compression);
//    }
//    this->mwTransformUp();
//    this->cropTree();
//    this->squareNorm = this->calcSquareNorm();
}

template<int D>
void FunctionTree<D>::power(double d) {
    NOT_IMPLEMENTED_ABORT;
//    this->purgeGenNodes();
//    for (int i = 0; i < this->endNodeTable.size(); i++) {
//        MWNode<D> &node = *this->endNodeTable[i];
//        node.mwTransform(Reconstruction);
//        node.cvTransform(MWNode<D>::Forward);
//        node.getCoefs() = node.getCoefs().array().pow(d);
//        node.cvTransform(MWNode<D>::Backward);
//        node.mwTransform(Compression);
//    }
//    this->mwTransformUp();
//    this->cropTree();
//    this->squareNorm = this->calcSquareNorm();
}

template<int D>
void FunctionTree<D>::normalize() {
    double norm = sqrt(this->getSquareNorm());
    *this *= (1.0/norm);
}

template<int D>
void FunctionTree<D>::map(const RepresentableFunction<1> &func) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
FunctionTree<D>& FunctionTree<D>::operator*=(double c) {
    for (int i = 0; i < this->endNodeTable.size(); i++) {
        MWNode<D> &node = *this->endNodeTable[i];
        if (not node.hasCoefs()) MSG_FATAL("No coefs");
        double *coefs = node.getCoefs();
        for (int j = 0; j < node.getNCoefs(); j++) {
            coefs[j] *= c;
        }
        node.calcNorms();
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
    return *this;
}

template<int D>
FunctionTree<D>& FunctionTree<D>::operator *=(const FunctionTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
FunctionTree<D>& FunctionTree<D>::operator +=(const FunctionTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
FunctionTree<D>& FunctionTree<D>::operator -=(const FunctionTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void FunctionTree<D>::getEndValues(VectorXd &data) {
    int nNodes = this->getNEndNodes();
    int nCoefs = this->getTDim()*this->getKp1_d();
    data = VectorXd::Zero(nNodes*nCoefs);
    for (int n = 0; n < nNodes; n++) {
        MWNode<D> &node = getEndFuncNode(n);
        node.mwTransform(Reconstruction);
        node.cvTransform(Forward);
        const double *c = node.getCoefs();
        for (int i = 0; i < nCoefs; i++) {
            data(n*nCoefs + i) = c[i];
        }
        node.cvTransform(Backward);
        node.mwTransform(Compression);
    }
}

template<int D>
void FunctionTree<D>::setEndValues(VectorXd &data) {
    int nNodes = this->getNEndNodes();
    int nCoefs = this->getTDim()*this->getKp1_d();
    for (int i = 0; i < nNodes; i++) {
        MWNode<D> &node = getEndFuncNode(i);
        const double *c = data.segment(i*nCoefs, nCoefs).data();
        node.setCoefBlock(0, nCoefs, c);
        node.cvTransform(Backward);
        node.mwTransform(Compression);
        node.setHasCoefs();
        node.calcNorms();
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

template class FunctionTree<1>;
template class FunctionTree<2>;
template class FunctionTree<3>;
