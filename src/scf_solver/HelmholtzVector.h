#pragma once

#include "qmoperators.h"
#include "qmfunctions.h"

namespace mrchem {

class HelmholtzVector final {
public:
    HelmholtzVector(double build, double thrs = -1.0);
    ~HelmholtzVector() { }

    void setup(double prec, const DoubleVector &energies);
    void clear();

    void setThreshold(double thrs) { this->threshold = thrs; }
    double getThreshold() const { return this->threshold; }

    double getLambda(int i) const { return this->lambda[i]; }
    DoubleVector getLambda() const;

    mrcpp::HelmholtzOperator& operator[](int i);
    const mrcpp::HelmholtzOperator& operator[](int i) const;

    int printTreeSizes() const;

    Orbital operator()(int i, Orbital inp);
    OrbitalVector operator()(OrbitalVector &inp);

private:
    double threshold; //For re-using operators. Negative means always recreate
    double build_prec;
    double apply_prec;

    std::vector<int> oper_idx;
    std::vector<double> lambda;
    std::vector<mrcpp::HelmholtzOperator *> operators;

    int initHelmholtzOperator(double energy, int i);
    void clearUnused();
};

} //namespace mrchem
