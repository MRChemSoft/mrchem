#include "ConfinementFunction.h"
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

    return f;
  }

} // namespace mrchem
