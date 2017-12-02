#pragma once

#include "Polynomial.h"

class LegendrePoly: public Polynomial {
public:
    LegendrePoly(int k, double n = 1.0, double l = 0.0);
    virtual ~LegendrePoly() { }

    Eigen::Vector2d firstDerivative(double x) const;
    Eigen::Vector3d secondDerivative(double x) const;

protected:
    void computeLegendrePolynomial(int k);
};

