#pragma once

#include "QMPotential.h"
#include "Density.h"

namespace mrchem {

class CoulombPotential final : public QMPotential {
public:
    CoulombPotential(mrcpp::PoissonOperator &P, OrbitalVector &Phi);
    ~CoulombPotential() { }

    void setup(double prec);
    void clear();

protected:
    Density density;                    // Density that defines the potential
    OrbitalVector *orbitals;            // Pointer to external object
    mrcpp::PoissonOperator *poisson;    // Pointer to external object

    void setupDensity(double prec);
    void setupPotential(double prec);
};

} //namespace mrchem
