#ifndef QMFUNCTION_H
#define QMFUNCTION_H

#include <complex>

#include "constants.h"
#include "FunctionTree.h"

template<int D>
class QMFunction {
public:
    QMFunction(FunctionTree<D> *r = 0, FunctionTree<D> *i = 0) : re(r), im(i) { }
    QMFunction(const QMFunction &func) { NOT_IMPLEMENTED_ABORT; }
    QMFunction<D> &operator=(const QMFunction<D> &func);
    virtual ~QMFunction() { clearReal(true); clearImag(true); }

    bool hasReal() const { if (this->re == 0) return false; return true; }
    bool hasImag() const { if (this->im == 0) return false; return true; }

    void allocReal();
    void allocImag();

    void clearReal(bool free);
    void clearImag(bool free);

    void setReal(FunctionTree<D> *real) { this->re = real; }
    void setImag(FunctionTree<D> *imag) { this->im = imag; }

    FunctionTree<D> &real() { return *this->re; }
    FunctionTree<D> &imag() { return *this->im; }

    const FunctionTree<D> &real() const { return *this->re; }
    const FunctionTree<D> &imag() const { return *this->im; }

    int getNNodes(int type = Total) const;
    double getSquareNorm(int type = Total) const;
    std::complex<double> dot(QMFunction<D> &ket);

    void normalize();
    void operator*=(double c);

    enum NumType { Total, Real, Imag };

protected:
    FunctionTree<D> *re;
    FunctionTree<D> *im;
};

#endif // QMFUNCTION_H
