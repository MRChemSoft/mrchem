#pragma once

#include "GroundStateSolver.h"

class Accelerator;

class OrbitalOptimizer : public GroundStateSolver {
public:
    OrbitalOptimizer(HelmholtzOperatorSet &h,
                     Accelerator *k = 0);
    virtual ~OrbitalOptimizer();

    void setup(FockOperator &fock, OrbitalVector &phi, Eigen::MatrixXd &F);
    void clear();

    virtual bool optimize();

protected:
    Accelerator *kain; // Pointer to external object, do not delete!

    void printTreeSizes() const;
};

