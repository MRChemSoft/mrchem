#pragma once

#include <string>

#include "MWProjector.h"
#include "GridGenerator.h"

class OrbitalVector;
class OrbitalExp;
class Nuclei;

class OrbitalProjector {
public:
    OrbitalProjector(double prec, int max_scale);
    virtual ~OrbitalProjector() { }

    void setPrecision(double prec) { this->project.setPrecision(prec); }

    OrbitalVector* operator()(const Nuclei &nucs);

    void operator()(OrbitalVector &phi,
                    const std::string &bf,
                    const std::string &mo);
    void operator()(OrbitalVector &phi,
                    const std::string &bf,
                    const std::string &mo_a,
                    const std::string &mo_b);
protected:
    MWProjector<3> project;
    GridGenerator<3> grid;

    OrbitalExp* readOrbitalExpansion(const std::string &bf,
                                     const std::string &mo);
};

