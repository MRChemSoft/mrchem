#include "QMFunction.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;

template<int D>
QMFunction<D>& QMFunction<D>::operator=(const QMFunction<D> &func) {
    if (this != &func) {
        if (this->hasReal()) MSG_ERROR("Function not empty");
        if (this->hasImag()) MSG_ERROR("Function not empty");
        this->re = func.re;
        this->im = func.im;
    }
    return *this;
}

template<int D>
void QMFunction<D>::allocReal() {
    if (this->hasReal()) MSG_ERROR("Function not empty");
    this->re = new FunctionTree<D>(*MRA);
}

template<int D>
void QMFunction<D>::allocImag() {
    if (this->hasImag()) MSG_ERROR("Function not empty");
    this->im = new FunctionTree<D>(*MRA);
}

template<int D>
void QMFunction<D>::clearReal(bool free) {
    if (this->hasReal() and free) delete this->re;
    this->re = 0;
}

template<int D>
void QMFunction<D>::clearImag(bool free) {
    if (this->hasImag() and free) delete this->im;
    this->im = 0;
}

template<int D>
int QMFunction<D>::getNNodes(int type) const {
    int rNodes = 0;
    int iNodes = 0;
    if (this->hasReal()) rNodes = this->real().getNNodes();
    if (this->hasImag()) iNodes = this->imag().getNNodes();
    if (type == Real) return rNodes;
    if (type == Imag) return iNodes;
    return rNodes + iNodes;
}

template<int D>
double QMFunction<D>::getSquareNorm(int type) const {
    double rNorm = 0;
    double iNorm = 0;
    if (this->hasReal()) rNorm = this->real().getSquareNorm();
    if (this->hasImag()) iNorm = this->imag().getSquareNorm();
    if (type == Real) return rNorm;
    if (type == Imag) return iNorm;
    return rNorm + iNorm;
}

template<int D>
complex<double> QMFunction<D>::dot(QMFunction<D> &ket) {
    QMFunction<D> &bra = *this;
    double rDot = 0.0;
    double iDot = 0.0;
    if (bra.hasReal() and ket.hasReal()) rDot += bra.real().dot(ket.real());
    if (bra.hasImag() and ket.hasImag()) rDot += bra.imag().dot(ket.imag());
    if (bra.hasReal() and ket.hasImag()) iDot += bra.real().dot(ket.imag());
    if (bra.hasImag() and ket.hasReal()) iDot -= bra.imag().dot(ket.real());
    return complex<double>(rDot, iDot);
}

template<int D>
void QMFunction<D>::normalize() {
    double norm = sqrt(this->getSquareNorm(Total));
    *this *= 1.0/norm;
}

template<int D>
void QMFunction<D>::operator*=(double c) {
    if (hasReal()) this->real() *= c;
    if (hasImag()) this->imag() *= c;
}

template class QMFunction<3>;
