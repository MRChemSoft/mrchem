#pragma once

//#include <MRCPP/mrcpp_declarations.h>

#include "qmoperators/QMPotential.h"
#include <vector>
#include "chemistry/Molecule.h"

namespace mrchem {

  class ConfinementPotential : public QMPotential {

  public:
    ConfinementPotential(double r_0, const int N, double s, std::vector<double> R, std::vector<mrcpp::Coord<3>> centers);

    auto getRadius() { return this->radius; }
    auto getStiffness() { return this->stiffness; }
    auto getSlope() { return this->slope; }
    auto getCavityradii() { return this->cavity_radii; }
    auto getCenters() { return this->centers; }

  protected:

    const int stiffness; // Stiffness parameter
    double radius;
    double slope; // Slope of the transition betwee boundaries
    std::vector<double> cavity_radii;
    std::vector<mrcpp::Coord<3>> centers;

    void setup(double prec);
    void clear();

  };

} // namespace mrchem
