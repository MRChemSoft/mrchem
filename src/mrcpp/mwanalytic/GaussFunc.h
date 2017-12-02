/*
 *
 * \breif Tools to handle gaussian functions and/or expansions in gaussian functions.
 * Implemented as a template class in the dimensionality:
 *
 *- Monodimensional gaussian (Gaussian<1>):
 *
 * \f$ g(x) = c e^{-\alpha (x-x^0)^2} \f$
 *
 * - Multidimensional gaussian (Gaussian<d>):
 *
 * \f$ G(x) = \prod_{i=1}^d g_i(x_i)
 * = \prod_{i=1}^d \c_i e^{-\alpha_i (x_i-x_i^0)^2} \f$
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "Gaussian.h"

template<int D> class GaussPoly;

template<int D>
class GaussFunc: public Gaussian<D> {
public:
    GaussFunc(double alpha = 0.0, double coef = 1.0, const double pos[D] = 0,
            const int pow[D] = 0) : Gaussian<D>(alpha, coef, pos, pow) {}
    GaussFunc(const GaussFunc<D> &gf) : Gaussian<D>(gf) {}
    Gaussian<D> *copy() const;
    ~GaussFunc() { }

    double calcCoulombEnergy(GaussFunc<D> &gf);
    double calcSquareNorm();

    double evalf(const double *r) const;
    double evalf(double r, int dim) const;

    static double calcOverlap(GaussFunc<D> &a, GaussFunc<D> &b);
    double calcOverlap(GaussFunc<D> &b);
    double calcOverlap(GaussPoly<D> &b);

    GaussPoly<D> differentiate(int dir);

    void multInPlace(const GaussFunc<D> &g);
    void operator*=(const GaussFunc<D> &gf) { multInPlace(gf); }
    GaussPoly<D> mult(const GaussFunc<D> &g);
    GaussFunc<D> mult(double d);
    GaussPoly<D> operator*(const GaussFunc<D> &g) { return this->mult(g); }
    GaussFunc<D> operator*(double d) { return this->mult(d); }

    void setPower(int d, int power) {
        this->power[d] = power;
        this->squareNorm = -1.0;
    }
    void setPower(const int power[D]) {
        for (int i = 0; i < D; i++) {
            this->power[i] = power[i];
        }
        this->squareNorm = -1.0;
    }
protected:
    static double ObaraSaika_ab(int power_a, int power_b, double pos_a,
            double pos_b, double expo_a, double expo_b);

};

