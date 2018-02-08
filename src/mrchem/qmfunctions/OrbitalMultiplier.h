#pragma once

class Orbital;

class OrbitalMultiplier {
public:
    OrbitalMultiplier(double prec, int max_scale);
    virtual ~OrbitalMultiplier() { }

    void setPrecision(double prec);

    void operator()(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b);
    void adjoint(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b);

protected:
    double multPrec;    
    void calcRealPart(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b);
    void calcImagPart(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b, bool adjoint);
};

