#pragma once

#include <string>

#include "MRCPP/Printer"
#include "parallel.h"

class OrbitalVector;
class OrbitalExp;
class Nuclei;

class OrbitalProjector {
public:
    OrbitalProjector(double prec);
    virtual ~OrbitalProjector() { }

    void setPrecision(double prec) { projPrec = prec; }

    OrbitalVector* operator()(const Nuclei &nucs);

    void operator()(OrbitalVector &phi,
                    const std::string &bf,
                    const std::string &mo);
    void operator()(OrbitalVector &phi,
                    const std::string &bf,
                    const std::string &mo_a,
                    const std::string &mo_b);
protected:
    double projPrec;
    OrbitalExp* readOrbitalExpansion(const std::string &bf,
                                     const std::string &mo);
};

