#include "FunctionNode.h"
#include "FunctionTree.h"
#include "MathUtils.h"
#include "QuadratureCache.h"
#include "blas.h"

using namespace std;
using namespace Eigen;

/** Function evaluation.
  * Evaluate all polynomials defined on the node. */
template<int D>
double FunctionNode<D>::evalf(const double *r) {
    if (not this->hasCoefs()) MSG_ERROR("Evaluating node without coefs");

    //we make many copies of the genChildren calls
    //the chance that 2 different nodes are treated in the same copy is small and therefore
    //the different critical sections will not interfer too much
    if (this->isLeafNode() or this->lockX!=0) {
      this->lockX=1;

      int Ncritical=16;
      switch (this->serialIx%Ncritical) {
      case 0:
#pragma omp critical(g0)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 1:
#pragma omp critical(g1)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 2:
#pragma omp critical(g2)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 3:
#pragma omp critical(g3)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 4:
#pragma omp critical(g4)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 5:
#pragma omp critical(g5)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 6:
#pragma omp critical(g6)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 7:
#pragma omp critical(g7)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 8:
#pragma omp critical(g8)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 9:
#pragma omp critical(g9)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 10:
#pragma omp critical(g10)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 11:
#pragma omp critical(g11)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 12:
#pragma omp critical(g12)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 13:
#pragma omp critical(g13)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 14:
#pragma omp critical(g14)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      case 15:
#pragma omp critical(g15)
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
        break;
      default:
#pragma omp critical
        if(this->isLeafNode()){
          this->genChildren();
          this->giveChildrenCoefs();
        }
      }
      this->lockX=0;
    }

    int cIdx = this->getChildIndex(r);
    assert(this->children[cIdx] != 0);
    return getFuncChild(cIdx).evalScaling(r);
}

template<int D>
double FunctionNode<D>::evalScaling(const double *r) const {
    if (not this->hasCoefs()) MSG_ERROR("Evaluating node without coefs");

    double arg[D];
    double n_factor = pow(2.0, this->getScale());
    const int *l_factor = this->getTranslation();
    for (int i = 0; i < D; i++) {
        arg[i] = r[i] * n_factor - (double) l_factor[i];
    }

    int fact[D + 1];
    for (int i = 0; i < D + 1; i++) {
        fact[i] = MathUtils::ipow(this->getKp1(), i);
    }

    MatrixXd val(this->getKp1(), D);
    const ScalingBasis &basis = this->getMWTree().getMRA().getScalingBasis();
    basis.evalf(arg, val);

    double result = 0.0;
//#pragma omp parallel for shared(fact) reduction(+:result)
    for (int i = 0; i < this->getKp1_d(); i++) {
        double temp = this->coefs[i];
        for (int j = 0; j < D; j++) {
            int k = (i % fact[j + 1]) / fact[j];
            temp *= val(k, j);
        }
        result += temp;
    }
    double n = (D * this->getScale()) / 2.0;
    double two_n = pow(2.0, n);
    return two_n * result;
}


/** Function integration.
  *
  * Wrapper for function integration, that requires different methods depending
  * on scaling type. Integrates the function represented on the node on the
  * full support of the node. */
template<int D>
double FunctionNode<D>::integrate() const {
    if (not this->hasCoefs()) {
        return 0.0;
    }
    switch (this->getScalingType()) {
    case Legendre:
        return integrateLegendre();
        break;
    case Interpol:
        return integrateInterpolating();
        break;
    default:
        MSG_FATAL("Invalid scalingType");
    }
}

/** Function integration, Legendre basis.
  *
  * Integrates the function represented on the node on the full support of the
  * node. The Legendre basis is particularly easy to integrate, as the work is
  * already done when calculating its coefficients. The coefficients of the
  * node is defined as the projection integral
  *          s_i = int f(x)phi_i(x)dx
  * and since the first Legendre function is the constant 1, the first
  * coefficient is simply the integral of f(x). */
template<int D>
double FunctionNode<D>::integrateLegendre() const {
    double n = (D * this->getScale()) / 2.0;
    double two_n = pow(2.0, -n);
    return two_n * this->getCoefs()[0];
}

/** Function integration, Interpolating basis.
  *
  * Integrates the function represented on the node on the full support of the
  * node. A bit more involved than in the Legendre basis, as is requires some
  * coupling of quadrature weights. */
template<int D>
double FunctionNode<D>::integrateInterpolating() const {
    int qOrder = this->getKp1();
    getQuadratureCache(qc);
    const VectorXd &weights = qc.getWeights(qOrder);

    double sqWeights[qOrder];
    for (int i = 0; i < qOrder; i++) {
        sqWeights[i] = sqrt(weights[i]);
    }

    int kp1_p[D];
    for (int i = 0; i < D; i++) {
        kp1_p[i] = MathUtils::ipow(qOrder, i);
    }

    VectorXd coefs;
    this->getCoefs(coefs);
    for (int p = 0; p < D; p++) {

        int n = 0;
        for (int i = 0; i < kp1_p[D - p - 1]; i++) {
            for (int j = 0; j < qOrder; j++) {
                for (int k = 0; k < kp1_p[p]; k++) {
                    coefs(n) *= sqWeights[j];
                    n++;
                }
            }
        }
    }
    double n = (D * this->getScale()) / 2.0;
    double two_n = pow(2.0, -n);
    double sum = coefs.segment(0, this->getKp1_d()).sum();

    return two_n * sum;
}

/** Inner product of the functions represented by the scaling basis of the nodes.
  *
  * Integrates the product of the functions represented by the scaling basis on
  * the node on the full support of the nodes. The scaling basis is fully
  * orthonormal, and the inner product is simply the dot product of the
  * coefficient vectors. Assumes the nodes have identical support. */
template<int D>
double FunctionNode<D>::dotScaling(const FunctionNode<D> &ket) const {
    const FunctionNode<D> &bra = *this;

    assert(bra.hasCoefs());
    assert(ket.hasCoefs());

    const double *a = bra.getCoefs();
    const double *b = ket.getCoefs();

    int size = bra.getKp1_d();
#ifdef HAVE_BLAS
    return ddot(size, a, 1, b, 1);
#else
    double result = 0.0;
    for (int i = 0; i < size; i++) {
        result += a[i]*b[i];
    }
    return result;
#endif
}

/** Inner product of the functions represented by the wavelet basis of the nodes.
  *
  * Integrates the product of the functions represented by the wavelet basis on
  * the node on the full support of the nodes. The wavelet basis is fully
  * orthonormal, and the inner product is simply the dot product of the
  * coefficient vectors. Assumes the nodes have identical support. */
template<int D>
double FunctionNode<D>::dotWavelet(const FunctionNode<D> &ket) const {
    const FunctionNode<D> &bra = *this;
    if (bra.isGenNode() or ket.isGenNode()) {
        return 0.0;
    }

    assert(bra.hasCoefs());
    assert(ket.hasCoefs());

    const double *a = bra.getCoefs();
    const double *b = ket.getCoefs();

    int start = bra.getKp1_d();
    int size = (bra.getTDim() - 1) * start;
#ifdef HAVE_BLAS
    return ddot(size, &a[start], 1, &b[start], 1);
#else
    double result = 0.0;
    for (int i = 0; i < size; i++) {
        result += a[start+i]*b[start+i];
    }
    return result;
#endif
}

template<int D>
void FunctionNode<D>::setValues(const VectorXd &vec) {
    this->zeroCoefs();
    this->setCoefBlock(0, vec.size(), vec.data());
    this->cvTransform(Backward);
    this->mwTransform(Compression);
    this->setHasCoefs();
    this->calcNorms();
}

template<int D>
void FunctionNode<D>::getValues(VectorXd &vec) {
    vec = VectorXd::Zero(this->n_coefs);
    this->mwTransform(Reconstruction);
    this->cvTransform(Forward);
    for (int i = 0; i < this->n_coefs; i++) {
        vec(i) = this->coefs[i];
    }
    this->cvTransform(Backward);
    this->mwTransform(Compression);
}

template class FunctionNode<1>;
template class FunctionNode<2>;
template class FunctionNode<3>;
