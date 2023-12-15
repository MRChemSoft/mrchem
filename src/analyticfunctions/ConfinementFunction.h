#pragma once

#include <MRCPP/MWFunctions>

namespace mrchem {

  class ConfinementFunction : public mrcpp::RepresentableFunction<3> {
  public:
    ConfinementFunction(double r_0, const int N, double slope, std::vector<double> R, std::vector<mrcpp::Coord<3>> centers);

    auto getRadius() { return this->radius; }
    auto getStiffness() { return this->stiffness; }
    auto getSlope() { return this->slope; }
    auto getCavityradii() { return this->cavity_radii; }
    auto getCenters() { return this->centers; }
    
  protected:

    double radius; //
    const int stiffness; // stffness parameter
    double slope; // cavity boundary slope
    std::vector<double> cavity_radii;
    std::vector<mrcpp::Coord<3>> centers;

    double evalf(const mrcpp::Coord<3> &r) const override;

  };

} // namespace mrchem
