#pragma once

#include "GroundStateSolver.h"

namespace mrchem {

class EnergyOptimizer final : public GroundStateSolver {
public:
    EnergyOptimizer(HelmholtzVector &h);
    ~EnergyOptimizer();

    void setup(FockOperator &fock, OrbitalVector &phi, ComplexMatrix &F,
               FockOperator &fock_np1, OrbitalVector &phi_np1);
    void clear();

    bool optimize();

protected:
    FockOperator *fOper_np1;
    OrbitalVector *orbitals_np1;

    ComplexMatrix calcFockMatrixUpdate(double prec, OrbitalVector &dPhi_n);
};

} //namespace mrchem
