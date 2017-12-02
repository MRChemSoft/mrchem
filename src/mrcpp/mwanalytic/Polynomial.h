/**
 *
 *  Base class for general polynomials with reasonably advanced
 * properties. The Polynomial class(es) are not implemented in the
 * most efficient manner, because they are only evaluated a fixed
 * number of times in a few predefined points, and all other
 * evaluations are done by linear transformations. PolynomialCache
 * implements the fast, and static const versions of the various
 * 4Polynomials.
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "RepresentableFunction.h"

class Polynomial: public RepresentableFunction<1> {
public:
    Polynomial(int k = 0, const double *a = 0, const double *b = 0);
    Polynomial(const Eigen::VectorXd &c, const double *a=0, const double *b=0);
    Polynomial(double c, int k = 0, const double *a=0, const double *b=0);
    Polynomial(const Polynomial &poly);
    Polynomial &operator=(const Polynomial &poly);
    virtual ~Polynomial() { }

    double evalf(double x) const;
    double evalf(const double *r) const { return evalf(r[0]); }

    double getScaledLowerBound() const;
    double getScaledUpperBound() const;

    void normalize();
    double calcSquareNorm();

    double getTranslation() const { return this->L; }
    double getDilation() const { return this->N; }

    void setDilation(double n) { this->N = n; }
    void setTranslation(double l) { this->L = l; }
    void rescale(double n, double l);

    int size() const { return this->coefs.size(); } ///< Length of coefs vector
    int getOrder() const;
    void clearCoefs() { this->coefs = Eigen::VectorXd::Zero(1); }
    void setZero() { this->coefs = Eigen::VectorXd::Zero(this->coefs.size()); }
    void setCoefs(const Eigen::VectorXd &c) { this->coefs = c; }

    Eigen::VectorXd &getCoefs() { return this->coefs; }
    const Eigen::VectorXd &getCoefs() const { return this->coefs; }

    Polynomial calcDerivative() const;
    Polynomial calcAntiDerivative() const;

    void calcDerivativeInPlace();
    void calcAntiDerivativeInPlace();

    double integrate(const double *a = 0, const double *b = 0) const;
    double innerProduct(const Polynomial &p) const;

    void addInPlace(double c, const Polynomial &Q);
    Polynomial add(double c, const Polynomial &Q) const;

    Polynomial operator*(double c) const;
    Polynomial operator*(const Polynomial &Q) const;
    Polynomial operator+(const Polynomial &Q) const { return add(1.0, Q); }
    Polynomial operator-(const Polynomial &Q) const { return add(-1.0, Q); }
    Polynomial &operator*=(double c);
    Polynomial &operator*=(const Polynomial &Q);
    Polynomial &operator+=(const Polynomial &Q);
    Polynomial &operator-=(const Polynomial &Q);
protected:
    double N; ///< Dilation coeff
    double L; ///< Translation coeff
    Eigen::VectorXd coefs; ///< Expansion coefficients
};
