#pragma once

#include "GroundStateSolver.h"

namespace mrchem {

class Accelerator;

class OrbitalOptimizer final : public GroundStateSolver {
public:
    OrbitalOptimizer(HelmholtzVector &h, Accelerator *k = 0);
    ~OrbitalOptimizer();

    void setup(FockOperator &fock, OrbitalVector &phi, Eigen::MatrixXd &F);
    void clear();

    bool optimize();

protected:
    Accelerator *kain; // Pointer to external object, do not delete!

    void printTreeSizes() const;
};

} //namespace mrchem
