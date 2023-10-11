#pragma once

//#include <MRCPP/mrcpp_declarations.h>

#include "qmoperators/QMPotential.h"
#include <vector>

namespace mrchem {

  class ConfinementPotential : public QMPotential {
  public:
    ConfinementPotential(const double r_0, const int N);

    auto getRadius() { return this->radius; }
    auto getParam() { return this->stiffness; }

  protected:

    const double radius; // Cavity radius
    const int stiffness; // Stiffness parameter

    void setup(double prec);
    void clear();

  };

} // namespace mrchem
