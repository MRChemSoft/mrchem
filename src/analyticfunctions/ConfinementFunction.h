#pragma once

#include <MRCPP/MWFunctions>

namespace mrchem {

  class ConfinementFunction : public mrcpp::RepresentableFunction<3> {
  public:
    ConfinementFunction(const double r_0, const int N);

    auto getRadius() { return this->radius; }
    auto getParam() { return this->stiffness; }
    
  protected:

    const double radius; //cavity radius
    const int stiffness; // stffness parameter
    
    double evalf(const mrcpp::Coord<3> &r) const override;

  };

} // namespace mrchem
