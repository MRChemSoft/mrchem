/*
 *  \date Feb 22, 2018
 *  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 *
 *  \breif Global objects and variables
 */

#pragma once

#include <complex>

#include <Eigen/Core>

#include "MRCPP/MWFunctions"
#include "config.h"

class Getkw;

// clang-format off
namespace mrchem {

namespace PHYSCONST {
const double kJ      = 2625.49962;         // AU -> kJ/mol conversion factor
const double kcal    =  627.509469;        // AU -> kcal/mol conversion factor
const double eV      =   27.21138505;      // AU -> eV conversion factor
const double JT_m2   =   78.9451185;       // AU -> 10^(-30) J/T^2 (SI magnetizability)
const double Debye   =    2.54174623105;   // AU -> Debye conversion factor
const double alpha   =    7.2973525664;    // Fine structure constant
const double g_e     =    2.0023193043618; // Free-electron g-value
const double beta_e  =    0.5000000;       // Bohr magneton
const double beta_N  =    0.0002723;       // Nuclear magneton
} // namespace PHYSCONST

namespace MATHCONST {
const double pi      =    3.1415926535897932384;
const double sqrt_pi =    1.7724538509055160273;
} // namespace MATHCONST

namespace SPIN { enum type { Paired, Alpha, Beta }; }
namespace NUMBER { enum type { Total, Real, Imag }; }
namespace DENSITY { enum type { Total, Spin, Alpha, Beta }; }
namespace MRDFT {enum type {Function, Gradient, Hessian}; }

 
using ComplexInt = std::complex<int>;
using ComplexDouble = std::complex<double>;

using IntVector = Eigen::VectorXi;
using DoubleVector = Eigen::VectorXd;
using ComplexVector = Eigen::VectorXcd;

using IntMatrix = Eigen::MatrixXi;
using DoubleMatrix = Eigen::MatrixXd;
using ComplexMatrix = Eigen::MatrixXcd;

extern Getkw Input;
extern mrcpp::MultiResolutionAnalysis<3> *MRA; //< Global MRA

} //namespace mrchem
// clang-format off
