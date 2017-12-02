#pragma once

#include "MWAdder.h"
#include "MWMultiplier.h"
#include "GridGenerator.h"

class OrbitalVector;
class Orbital;
class Density;

class DensityProjector {
public:
    DensityProjector(double prec, int max_scale);
    virtual ~DensityProjector() { }

    void setPrecision(double prec);

    void operator()(Density &rho, Orbital &phi);
    void operator()(Density &rho, OrbitalVector &phi);

protected:
    MWAdder<3> add;
    MWMultiplier<3> mult;
    GridGenerator<3> grid;
};

