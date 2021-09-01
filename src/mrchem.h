/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#pragma once

#include <complex>

#include <Eigen/Core>

#include "MRCPP/MWFunctions"
#include "config.h"
#include "mrenum.h"

// clang-format off
namespace mrchem {

namespace PHYSCONST {
const double kJ        = 2625.49962;         // AU -> kJ/mol conversion factor
const double kcal      =  627.509469;        // AU -> kcal/mol conversion factor
const double eV        =   27.21138505;      // AU -> eV conversion factor
const double JT_m2     =   78.9451185;       // AU -> 10^(-30) J/T^2 (SI magnetizability)
const double Debye     =    2.54174623105;   // AU -> Debye conversion factor
const double cm_m1     =    1.0e7 / 45.5640; // AU -> cm^(-1) frequency -> wavelenght
const double alpha     =    7.2973525664;    // Fine structure constant
const double g_e       =    2.0023193043618; // Free-electron g-value
const double beta_e    =    0.5000000;       // Bohr magneton
const double beta_N    =    0.0002723;       // Nuclear magneton
const double alpha_inv =  137.035999084;     // Inverse of the fine structure constant alpha
} // namespace PHYSCONST

namespace MATHCONST {
const double pi      =    3.1415926535897932384;
const double sqrt_pi =    1.7724538509055160273;
} // namespace MATHCONST

namespace SPIN { enum type { Paired, Alpha, Beta }; }
namespace NUMBER { enum type { Total, Real, Imag }; }
namespace MRDFT {enum type {Function, Gradient, Hessian}; }

using ComplexInt = std::complex<int>;
using ComplexDouble = std::complex<double>;

using IntVector = Eigen::VectorXi;
using DoubleVector = Eigen::VectorXd;
using ComplexVector = Eigen::VectorXcd;

using IntMatrix = Eigen::MatrixXi;
using DoubleMatrix = Eigen::MatrixXd;
using ComplexMatrix = Eigen::MatrixXcd;

extern mrcpp::MultiResolutionAnalysis<3> *MRA; //< Global MRA

} //namespace mrchem
// clang-format off
