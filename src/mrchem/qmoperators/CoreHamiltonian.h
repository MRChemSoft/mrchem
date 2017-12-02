#pragma once

#include "FockOperator.h"

class CoreHamiltonian : public FockOperator {
public:
    CoreHamiltonian(KineticOperator &t,
                    NuclearPotential &v)
        : FockOperator(&t, &v) { }
    virtual ~CoreHamiltonian() { }
};

