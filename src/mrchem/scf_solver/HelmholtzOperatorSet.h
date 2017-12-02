#pragma once

#pragma GCC system_header
#include <Eigen/Core>
#include <vector>

#include "HelmholtzOperator.h"
#include "MWConvolution.h"

class Orbital;
class OrbitalVector;

class HelmholtzOperatorSet {
public:
    HelmholtzOperatorSet(double build, double thrs = -1.0);
    virtual ~HelmholtzOperatorSet() { clear(); }

    void setup(double prec, const Eigen::VectorXd &energies);
    void clear();

    void setThreshold(double thrs) { this->threshold = thrs; }
    double getThreshold() const { return this->threshold; }

    double getLambda(int i) const { return this->lambda[i]; }
    Eigen::VectorXd getLambda() const;
    HelmholtzOperator &getOperator(int i);
    printTreeSizes() const;

    void operator()(int i, Orbital &out, Orbital &inp);
    void operator()(OrbitalVector &out, OrbitalVector &inp);

private:
    double threshold; //For re-using operators. Negative means always recreate
    double build_prec;
    MWConvolution<3> apply;

    std::vector<int> operIdx;
    std::vector<double> lambda;
    std::vector<HelmholtzOperator *> operators;

    int initHelmholtzOperator(double energy, int i);
    void clearUnused();
};


