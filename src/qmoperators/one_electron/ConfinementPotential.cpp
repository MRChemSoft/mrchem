#include "ConfinementPotential.h"
#include "analyticfunctions/ConfinementFunction.h"

#include "utils/math_utils.h"
#include "chemistry/Molecule.h"

namespace mrchem {

  ConfinementPotential::ConfinementPotential(double r_0, const int N, double s, std::vector<double> R, std::vector<mrcpp::Coord<3>> coords)
    : QMPotential(1, false)
    , radius(r_0)
    , stiffness(N)
    , slope(s)
    , cavity_radii(R)
    , centers(coords) {}


  void ConfinementPotential::setup(double prec) {
    // Initalize the function representing the confinement
    ConfinementFunction *f_loc = nullptr;

    f_loc = new ConfinementFunction(this->getRadius(), this->getStiffness(), this->getSlope(), this->getCavityradii(), this->getCenters());

    // Project the potential onto the function representation
    mrcpp::ComplexFunction V_loc(false);
    mrcpp::cplxfunc::project(V_loc, *f_loc, NUMBER::Real, prec);
    
    mrcpp::ComplexFunction &V_c = (*this);
    setApplyPrec(prec);
    V_c = V_loc;
    delete f_loc;

  }

  void ConfinementPotential::clear() {
    mrcpp::ComplexFunction::free(NUMBER::Total);
    clearApplyPrec();
  }


} // namespace mrchem
