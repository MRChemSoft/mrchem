#include "lebedev_utils.hpp"

#include <cmath>
#include <tuple>

#include <Eigen/Core>

namespace detail {
auto
a1(double w) -> std::tuple<Eigen::Matrix<double, 6, 1>, Eigen::Matrix<double, 3, 6>>
{
  auto a = 1.0;

  Eigen::Matrix<double, 3, 6> xs{ { +a, -a, 0.0, 0.0, 0.0, 0.0 },
                                  { 0.0, 0.0, +a, -a, 0.0, 0.0 },
                                  { 0.0, 0.0, 0.0, 0.0, +a, -a } };

  return { Eigen::Matrix<double, 6, 1>::Constant(w), xs };
}

auto
a2(double w) -> std::tuple<Eigen::Matrix<double, 12, 1>, Eigen::Matrix<double, 3, 12>>
{
  auto a = 1.0 / std::sqrt(2);

  Eigen::Matrix<double, 3, 12> xs{ { +a, +a, -a, -a, +a, +a, -a, -a, 0., 0., 0., 0. },
                                   { +a, -a, +a, -a, 0., 0., 0., 0., +a, +a, -a, -a },
                                   { 0., 0., 0., 0., +a, -a, +a, -a, +a, -a, +a, -a } };

  return { Eigen::Matrix<double, 12, 1>::Constant(w), xs };
}

auto
a3(double w) -> std::tuple<Eigen::Matrix<double, 8, 1>, Eigen::Matrix<double, 3, 8>>
{
  auto a = 1.0 / std::sqrt(3);

  Eigen::Matrix<double, 3, 8> xs{ { +a, +a, +a, +a, -a, -a, -a, -a },
                                  { +a, +a, -a, -a, +a, +a, -a, -a },
                                  { +a, -a, +a, -a, +a, -a, +a, -a } };

  return { Eigen::Matrix<double, 8, 1>::Constant(w), xs };
}

auto
pq0(double w, double x) -> std::tuple<Eigen::Matrix<double, 24, 1>, Eigen::Matrix<double, 3, 24>>
{
  auto a = std::sin(x * M_PI);
  auto b = std::cos(x * M_PI);

  Eigen::Matrix<double, 3, 24> xs{
    { +a, -a, -a, +a, +b, -b, -b, +b, +a, -a, -a, +a, +b, -b, -b, +b, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { +b, +b, -b, -b, +a, +a, -a, -a, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, +a, -a, -a, +a, +b, -b, -b, +b },
    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, +b, +b, -b, -b, +a, +a, -a, -a, +b, +b, -b, -b, +a, +a, -a, -a }
  };

  return { Eigen::Matrix<double, 24, 1>::Constant(w), xs };
}

auto
llm(double w, double x) -> std::tuple<Eigen::Matrix<double, 24, 1>, Eigen::Matrix<double, 3, 24>>
{
  auto a = std::sin(x * M_PI) / std::sqrt(2);
  auto b = std::cos(x * M_PI);

  Eigen::Matrix<double, 3, 24> xs{
    { +a, -a, +a, -a, +a, -a, +a, -a, +a, -a, +a, -a, +a, -a, +a, -a, +b, +b, +b, +b, -b, -b, -b, -b },
    { +a, +a, -a, -a, +a, +a, -a, -a, +b, +b, +b, +b, -b, -b, -b, -b, +a, -a, +a, -a, +a, -a, +a, -a },
    { +b, +b, +b, +b, -b, -b, -b, -b, +a, +a, -a, -a, +a, +a, -a, -a, +a, +a, -a, -a, +a, +a, -a, -a }
  };

  return { Eigen::Matrix<double, 24, 1>::Constant(w), xs };
}

auto
rsw(double w, double x, double y) -> std::tuple<Eigen::Matrix<double, 48, 1>, Eigen::Matrix<double, 3, 48>>
{
  auto sinx = std::sin(x * M_PI);
  auto siny = std::sin(y * M_PI);
  auto cosx = std::cos(x * M_PI);
  auto cosy = std::cos(y * M_PI);

  auto a = sinx * cosy;
  auto b = sinx * siny;
  auto c = cosx;

  Eigen::Matrix<double, 3, 48> xs{
    { +a, +c, +b, +b, +c, +a, -a, +c, +b, +b, +c, -a, +a, +c, -b, -b, +c, +a, +a, -c, +b, +b, -c, +a,
      -a, +c, -b, -b, +c, -a, -a, -c, +b, +b, -c, -a, +a, -c, -b, -b, -c, +a, -a, -c, -b, -b, -c, -a },
    { +b, +a, +c, +a, +b, +c, +b, -a, +c, -a, +b, +c, -b, +a, +c, +a, -b, +c, +b, +a, -c, +a, +b, -c,
      -b, -a, +c, -a, -b, +c, +b, -a, -c, -a, +b, -c, -b, +a, -c, +a, -b, -c, -b, -a, -c, -a, -b, -c },
    { +c, +b, +a, +c, +a, +b, +c, +b, -a, +c, -a, +b, +c, -b, +a, +c, +a, -b, -c, +b, +a, -c, +a, +b,
      +c, -b, -a, +c, -a, -b, -c, +b, -a, -c, -a, +b, -c, -b, +a, -c, +a, -b, -c, -b, -a, -c, -a, -b }
  };

  return { Eigen::Matrix<double, 48, 1>::Constant(w), xs };
}

auto
n_points_to_degree(size_t N) -> size_t
{
  switch (N) {
    case 6:
      return 3;
    case 14:
      return 5;
    case 26:
      return 7;
    case 38:
      return 9;
    case 50:
      return 11;
    case 74:
      return 13;
    case 86:
      return 15;
    case 110:
      return 17;
    case 146:
      return 19;
    case 170:
      return 21;
    case 194:
      return 23;
    case 230:
      return 25;
    case 266:
      return 27;
    case 302:
      return 39;
    case 350:
      return 31;
    case 434:
      return 35;
    case 590:
      return 41;
    case 770:
      return 47;
    case 974:
      return 53;
    case 1202:
      return 59;
    case 1454:
      return 65;
    case 1730:
      return 71;
    case 2030:
      return 77;
    case 2354:
      return 83;
    case 2702:
      return 89;
    case 3074:
      return 95;
    case 3470:
      return 101;
    case 3890:
      return 107;
    case 4334:
      return 113;
    case 4802:
      return 119;
    case 5294:
      return 125;
    case 5810:
      return 131;
    default:
      // FIXME error handling
      std::abort();
  }
}
} // namespace detail
