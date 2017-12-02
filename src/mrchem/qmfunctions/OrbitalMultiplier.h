#pragma once

#include "MWAdder.h"
#include "MWMultiplier.h"
#include "GridGenerator.h"

class Orbital;

class OrbitalMultiplier {
public:
    OrbitalMultiplier(double prec, int max_scale);
    virtual ~OrbitalMultiplier() { }

    void setPrecision(double prec);

    void operator()(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b);
    void adjoint(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b);

protected:
    MWAdder<3> add;
    MWMultiplier<3> mult;
    GridGenerator<3> grid;

    void calcRealPart(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b);
    void calcImagPart(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b, bool adjoint);
};

