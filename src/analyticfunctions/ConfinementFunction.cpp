#include "ConfinementFunction.h"
#include "environment/Cavity.h"
#include "utils/math_utils.h"

namespace mrchem {

  ConfinementFunction::ConfinementFunction(const double r_0, const int N)
    : radius(r_0)
    , stiffness(N) {}

  double ConfinementFunction::evalf(const mrcpp::Coord<3> &r) const {
    double f = 0.0;
    const mrcpp::Coord<3> center = {0.0, 0.0, 0.0};
    const int N = this->stiffness;
    auto r_0 = this->radius;

    f += std::pow(math_utils::calc_distance(r, center) / r_0, N);

    std::vector<mrcpp::Coord<3>> coords = {{0.0, 0.0, 0.0}};
    std::vector<double> R = {1.0};
    double slope = 0.2;

    Cavity sphere(coords, R, slope);

    f * (1 - sphere.evalf(r) );

    return f;
  }

} // namespace mrchem
