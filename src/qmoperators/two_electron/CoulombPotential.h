#pragma once

#include "CoulombOperator.h"

class CoulombPotential : public CoulombOperator {
public:
    CoulombPotential(PoissonOperator &P, OrbitalVector &phi)
        : CoulombOperator(P, phi) { }
    virtual ~CoulombPotential() { }

    virtual void setup(double prec) {
        setApplyPrec(prec);
        calcDensity(this->density, *this->orbitals);
        calcPotential(*this, this->density);
    }
    virtual void clear() {
        clearReal(true);
        clearImag(true);
        this->density.clear();
        clearApplyPrec();
    }

protected:
    void calcDensity(Density &rho, OrbitalVector &phi);
    void calcPotential(QMPotential &V, Density &rho);
};

