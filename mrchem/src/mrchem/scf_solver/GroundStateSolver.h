#pragma once

#include "SCF.h"
#include "properties/SCFEnergy.h"

/** @class GroundStateSolver
 *
 * @brief Abstract class for different types of ground state SCF solvers
 *
 * The ground state SCF solvers share some common features which are collected in this
 * abstract base class. This is mainly the construction of the argument that is used
 * for the Helmholtz operators.
 */

namespace mrchem {

class GroundStateSolver : public SCF {
public:
    GroundStateSolver(HelmholtzVector &h);

protected:
    std::vector<SCFEnergy> energy;

    ComplexMatrix *fMat_n;          ///< Fock matrix (pointer to external object)
    FockOperator  *fOper_n;         ///< Fock operator (pointer to external object)
    OrbitalVector *orbitals_n;      ///< Orbtials (pointer to external object)

    OrbitalVector setupHelmholtzArguments(FockOperator &fock,
                                          const ComplexMatrix &M,
                                          OrbitalVector &Phi,
                                          bool clearFock);
    void printProperty() const;
    double calcProperty();
    double calcPropertyError() const;
};

} //namespace mrchem
