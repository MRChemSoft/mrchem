#include "ConfinementFunction.h"
#include "environment/Cavity.h"
#include "utils/math_utils.h"

namespace mrchem {

  ConfinementFunction::ConfinementFunction(double r_0, const int N, double s, std::vector<double> R, std::vector<mrcpp::Coord<3>> centers)
    : radius(r_0)
    , slope(s)
    , stiffness(N)
    , cavity_radii(R)
    , centers(centers) {}

  double ConfinementFunction::evalf(const mrcpp::Coord<3> &r) const {
    double f = 0.0;

    const int N = this->stiffness;
    auto r_0 = this->radius;
    auto s = this->slope;
    auto R = this->cavity_radii;
    auto coords = this->centers;

    Cavity sphere(coords, R, s);

    //double min_dist = 1000.0;
    //double &rad = min_dist;
    //for (auto& center : coords) {
    //  double distance = math_utils::calc_distance(r, center);
    // if (distance < min_dist) {
    //rad = distance;
    //  }
    //}

    //f += std::pow(min_dist / r_0, N);
    mrcpp::Coord<3> origin = {0.0, 0.0, 0.0};
    f += std::pow(math_utils::calc_distance(origin, r) / r_0, N);
    f *= (1 - sphere.evalf(r));

    return f;
  }

} // namespace mrchem
