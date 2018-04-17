#pragma once

#include "QMPotential.h"
#include "Density.h"

/** @class CoulombPotential
 *
 * @brief Coulomb potential defined by a particular set of orbitals
 *
 * The OrbitalVector defining the operator is fixed throughout the operators life time,
 * but the orbitals themselves are allowed to change in between each application. When
 * the operator has been setup, it will be fixed until it's cleared and setup again
 * (even if the orbitals change).
 */

namespace mrchem {

class CoulombPotential final : public QMPotential {
public:
    CoulombPotential(mrcpp::PoissonOperator &P, OrbitalVector &Phi);
    ~CoulombPotential() { }

protected:
    Density density;                        // Density that defines the potential
    OrbitalVector *orbitals;          // Pointer to external object
    mrcpp::PoissonOperator *poisson;  // Pointer to external object

    void setup(double prec);
    void clear();

    void setupDensity(double prec);
    void setupPotential(double prec);
};

} //namespace mrchem
