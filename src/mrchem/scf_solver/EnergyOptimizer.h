#pragma once

#include "GroundStateSolver.h"

class EnergyOptimizer : public GroundStateSolver {
public:
    EnergyOptimizer(HelmholtzOperatorSet &h);
    virtual ~EnergyOptimizer();

    void setup(FockOperator &fock, OrbitalVector &phi, Eigen::MatrixXd &F,
               FockOperator &fock_np1, OrbitalVector &phi_np1);
    void clear();

    virtual bool optimize();

protected:
    FockOperator *fOper_np1;

    Eigen::MatrixXd calcFockMatrixUpdate();
};

