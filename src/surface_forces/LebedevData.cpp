#include "LebedevData.h"

#include <algorithm>
#include <cmath>
#include <tuple>

#include <Eigen/Core>

#include "detail/Lebedev_003_0006.hpp"
#include "detail/Lebedev_005_0014.hpp"
#include "detail/Lebedev_007_0026.hpp"
#include "detail/Lebedev_009_0038.hpp"
#include "detail/Lebedev_011_0050.hpp"
#include "detail/Lebedev_013_0074.hpp"
#include "detail/Lebedev_015_0086.hpp"
#include "detail/Lebedev_017_0110.hpp"
#include "detail/Lebedev_019_0146.hpp"
#include "detail/Lebedev_021_0170.hpp"
#include "detail/Lebedev_023_0194.hpp"
#include "detail/Lebedev_025_0230.hpp"
#include "detail/Lebedev_027_0266.hpp"
#include "detail/Lebedev_029_0302.hpp"
#include "detail/Lebedev_031_0350.hpp"
#include "detail/Lebedev_035_0434.hpp"
#include "detail/Lebedev_041_0590.hpp"
#include "detail/Lebedev_047_0770.hpp"
#include "detail/Lebedev_053_0974.hpp"
#include "detail/Lebedev_059_1202.hpp"
#include "detail/Lebedev_065_1454.hpp"
#include "detail/Lebedev_071_1730.hpp"
#include "detail/Lebedev_077_2030.hpp"
#include "detail/Lebedev_083_2354.hpp"
#include "detail/Lebedev_089_2702.hpp"
#include "detail/Lebedev_095_3074.hpp"
#include "detail/Lebedev_101_3470.hpp"
#include "detail/Lebedev_107_3890.hpp"
#include "detail/Lebedev_113_4334.hpp"
#include "detail/Lebedev_119_4802.hpp"
#include "detail/Lebedev_125_5294.hpp"
#include "detail/Lebedev_131_5810.hpp"

auto
lebedev(size_t N) -> std::tuple<Eigen::VectorXd, Eigen::Matrix3Xd>
{
  std::array available = {
    6,   14,  26,  38,   50,   74,   86,   110,  146,  170,  194,  230,  266,  302,  350,  434,
    590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810
  };

  const auto lb = std::lower_bound(available.cbegin(), available.cend(), N);
  auto closest = (lb == available.cend()) ? available.back() : *lb;

  if (lb != available.cbegin()) {
    auto prec = lb - 1;
    if (std::abs<int>(closest - N) > std::abs<int>(*prec - N))
      closest = *prec;
  }

  switch (closest) {
    case 6:
      return detail::Lebedev_003_0006::quadrature();
    case 14:
      return detail::Lebedev_005_0014::quadrature();
    case 26:
      return detail::Lebedev_007_0026::quadrature();
    case 38:
      return detail::Lebedev_009_0038::quadrature();
    case 50:
      return detail::Lebedev_011_0050::quadrature();
    case 74:
      return detail::Lebedev_013_0074::quadrature();
    case 86:
      return detail::Lebedev_015_0086::quadrature();
    case 110:
      return detail::Lebedev_017_0110::quadrature();
    case 146:
      return detail::Lebedev_019_0146::quadrature();
    case 170:
      return detail::Lebedev_021_0170::quadrature();
    case 194:
      return detail::Lebedev_023_0194::quadrature();
    case 230:
      return detail::Lebedev_025_0230::quadrature();
    case 266:
      return detail::Lebedev_027_0266::quadrature();
    case 302:
      return detail::Lebedev_029_0302::quadrature();
    case 350:
      return detail::Lebedev_031_0350::quadrature();
    case 434:
      return detail::Lebedev_035_0434::quadrature();
    case 590:
      return detail::Lebedev_041_0590::quadrature();
    case 770:
      return detail::Lebedev_047_0770::quadrature();
    case 974:
      return detail::Lebedev_053_0974::quadrature();
    case 1202:
      return detail::Lebedev_059_1202::quadrature();
    case 1454:
      return detail::Lebedev_065_1454::quadrature();
    case 1730:
      return detail::Lebedev_071_1730::quadrature();
    case 2030:
      return detail::Lebedev_077_2030::quadrature();
    case 2354:
      return detail::Lebedev_083_2354::quadrature();
    case 2702:
      return detail::Lebedev_089_2702::quadrature();
    case 3074:
      return detail::Lebedev_095_3074::quadrature();
    case 3470:
      return detail::Lebedev_101_3470::quadrature();
    case 3890:
      return detail::Lebedev_107_3890::quadrature();
    case 4334:
      return detail::Lebedev_113_4334::quadrature();
    case 4802:
      return detail::Lebedev_119_4802::quadrature();
    case 5294:
      return detail::Lebedev_125_5294::quadrature();
    case 5810:
      return detail::Lebedev_131_5810::quadrature();

    default:
      // FIXME error handling
      std::abort();
  }
}