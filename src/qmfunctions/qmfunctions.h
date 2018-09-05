#pragma once

#include "mrchem.h"

namespace mrchem {

namespace SPIN { enum type { Paired, Alpha, Beta }; }
namespace NUMBER { enum type { Total, Real, Imag }; }
namespace DENSITY { enum type { Total, Spin, Alpha, Beta }; }

class QMFunction;
using QMFunctionVector = std::vector<std::tuple<double, QMFunction>>;

class Orbital;
using OrbitalChunk = std::vector<std::tuple<int, Orbital>>;
using OrbitalVector = std::vector<Orbital>;

class Density;
using DensityChunk = std::vector<std::tuple<int, Density>>;
using DensityVector = std::vector<Density>;

} //namespace mrchem
