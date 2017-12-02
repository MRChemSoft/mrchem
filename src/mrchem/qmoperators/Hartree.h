#pragma once

#include "FockOperator.h"

class Hartree : public FockOperator {
public:
    Hartree(KineticOperator &t,
            NuclearPotential &v,
            CoulombOperator &j)
        : FockOperator(&t, &v, &j) { }
    virtual ~Hartree() { }
};

