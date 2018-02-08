#pragma once

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
    double densPrec;
};

