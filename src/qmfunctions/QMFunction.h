#pragma once

#include "qmfunctions.h"

namespace mrchem {

class QMFunction {
public:
    QMFunction(mrcpp::FunctionTree<3> *r = nullptr,
               mrcpp::FunctionTree<3> *i = nullptr)
        : re(r), im(i) { }
    QMFunction(const QMFunction &func) : re(func.re), im(func.im) { }
    QMFunction &operator=(const QMFunction &func);
    virtual ~QMFunction() = default;

    void alloc(int type = NUMBER::Total);
    void clear(int type = NUMBER::Total);
    void free(int type = NUMBER::Total);

    int getNNodes(int type = NUMBER::Total) const;

    bool hasReal() const { return (this->re == nullptr) ? false : true; }
    bool hasImag() const { return (this->im == nullptr) ? false : true; }

    mrcpp::FunctionTree<3> &real() { return *this->re; }
    mrcpp::FunctionTree<3> &imag() { return *this->im; }

    const mrcpp::FunctionTree<3> &real() const { return *this->re; }
    const mrcpp::FunctionTree<3> &imag() const { return *this->im; }

    void setReal(mrcpp::FunctionTree<3> *real) { this->re = real; }
    void setImag(mrcpp::FunctionTree<3> *imag) { this->im = imag; }

protected:
    mrcpp::FunctionTree<3> *re;     ///< Real part of function
    mrcpp::FunctionTree<3> *im;     ///< Imaginary part of function
};

} // namespace mrchem
