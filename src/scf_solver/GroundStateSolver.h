#pragma once

#include "SCF.h"
#include "SCFEnergy.h"

namespace mrchem {

class GroundStateSolver : public SCF {
public:
    GroundStateSolver(HelmholtzVector &h);
    virtual ~GroundStateSolver();

protected:
    std::vector<SCFEnergy> energy;

    ComplexMatrix *fMat_n;
    FockOperator  *fOper_n;
    OrbitalVector *orbitals_n;

    OrbitalVector setupHelmholtzArguments(FockOperator &fock,
                                          const ComplexMatrix &M,
                                          OrbitalVector &Phi,
                                          bool adjoint = false,
                                          bool clearFock = false);
    void printProperty() const;
    double calcProperty();
    double calcPropertyError() const;
};

} //namespace mrchem
