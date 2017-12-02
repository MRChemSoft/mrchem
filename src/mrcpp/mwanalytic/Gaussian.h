/**
 *
 *  Base class for Gaussian type functions
 */

#pragma once

#include <iostream>
#include <cmath>

#pragma GCC system_header
#include <Eigen/Core>

#include "RepresentableFunction.h"

template<int D> class GaussFunc;
template<int D> class GaussPoly;
template<int D> class GaussExp;

template<int D>
class Gaussian: public RepresentableFunction<D> {
public:
    Gaussian(double a, double c, const double r[D], const int p[D]);
    virtual Gaussian<D> *copy() const = 0;
    virtual ~Gaussian();

    virtual double evalf(const double *r) const = 0;
    virtual double evalf(double r, int dim) const = 0;
    void evalf(const Eigen::MatrixXd &points, Eigen::MatrixXd &values) const;

    virtual double calcSquareNorm() = 0;
    virtual double calcOverlap(GaussFunc<D> &b) = 0;
    virtual double calcOverlap(GaussPoly<D> &b) = 0;

    virtual GaussPoly<D> differentiate(int dir) = 0;

    void calcScreening(double stdDeviations);

    double getSquareNorm() {
        if (this->squareNorm < 0.0) {
            calcSquareNorm();
        }
        return this->squareNorm;
    }
    void normalize() {
        double norm = sqrt(getSquareNorm());
        multConstInPlace(1.0/norm);
    }
    void multPureGauss(const Gaussian<D> &lhs, const Gaussian<D> &rhs);
    void multConstInPlace(double c) {
        this->coef *= c;
        this->calcSquareNorm();
    }
    void operator*=(double c) {
        multConstInPlace(c);
    }

    bool checkScreen(int n, const int *l) const;
    int getPower(int i) const { return power[i]; }
    bool getScreen() const { return screen; }
    const int *getPower() const { return power; }
    const double *getPos() const { return pos; }
    double getCoef() const { return coef; }
    double getExp() const { return alpha; }

    virtual void setPower(const int power[D]) = 0;
    virtual void setPower(int d, int power) = 0;
    void setScreen(bool _screen) { this->screen = _screen; }
    void setCoef(double cf) { this->coef = cf; this->squareNorm = -1.0; }
    void setExp(double _alpha) { this->alpha = _alpha; this->squareNorm = -1.0; }
    void setPos(const double r[D]) {
        for (int d = 0; d < D; d++) {
            this->pos[d] = r[d];
        }
        this->squareNorm = -1.0;
    }
    friend std::ostream& operator<<(std::ostream &o, const Gaussian<D> &gauss)
    {
        o << "Exp:   " << gauss.getExp() << std::endl;
        o << "Coef:  "<< gauss.getCoef() << std::endl;
        o << "Pos:   ";
        for (int i = 0; i < D; i++) {
            o << gauss.getPos()[i] << " ";
        }
        o << std::endl;
        o << "Power: ";
        for (int i = 0; i < D; i++) {
            o << gauss.getPower(i) << " ";
        }
        return o;
    }

    friend class GaussExp<D>;
protected:
    bool screen;
    int power[D];		/**< max power in each dim  */
    double coef;		/**< constant factor */
    double alpha;		/**< exponent  */
    double pos[D];		/**< center  */
    double squareNorm;

    bool isVisibleAtScale(int scale, int nQuadPts) const;
    bool isZeroOnInterval(const double *a, const double *b) const;
};

