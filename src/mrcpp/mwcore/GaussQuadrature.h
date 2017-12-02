#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "RepresentableFunction.h"

const int MaxGaussOrder = 42;
static const double EPS = 3.0e-12;
static const int NewtonMaxIter = 10;
static const int MaxQuadratureDim = 7;

class GaussQuadrature {
public:
    GaussQuadrature(int k, double a = -1.0, double b = 1.0, int inter = 1);
    virtual ~GaussQuadrature() { }

    double integrate(RepresentableFunction<1> &func) const;
    double integrate(RepresentableFunction<2> &func) const;
    double integrate(RepresentableFunction<3> &func) const;

    void setIntervals(int i);
    void setBounds(double a, double b);

    int getIntervals() const { return this->intervals; }
    double getUpperBound() const { return this->B; }
    double getLowerBound() const { return this->A; }

    const Eigen::VectorXd &getRoots() const { return this->roots; }
    const Eigen::VectorXd &getWeights() const { return this->weights; }
    const Eigen::VectorXd &getUnscaledRoots() const { return this->unscaledRoots; }
    const Eigen::VectorXd &getUnscaledWeights() const { return this->unscaledWeights; }

protected:
    int order;
    double A;
    double B;
    int intervals;
    int npts;
    Eigen::VectorXd roots;
    Eigen::VectorXd weights;
    Eigen::VectorXd unscaledRoots;
    Eigen::VectorXd unscaledWeights;

    void rescaleRoots(Eigen::VectorXd &rts, double a, double b, int inter = 1) const;
    void rescaleWeights(Eigen::VectorXd &wgts, double a, double b, int inter = 1) const;

    void calcScaledPtsWgts();
    int calcGaussPtsWgts();

    double integrate_nd(RepresentableFunction<3> &func, int axis = 0) const;
};

