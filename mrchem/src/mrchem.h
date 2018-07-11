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

#include "config.h"
#include "MRCPP/MWFunctions"

class Getkw;

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
}

namespace MATHCONST {
    const double pi      =    3.1415926535897932384;
    const double sqrt_pi =    1.7724538509055160273;
}

typedef std::complex<double> ComplexDouble;
typedef std::function<double (const double *r)> DoubleFunction;
typedef std::function<ComplexDouble (const double *r)> ComplexFunction;

typedef Eigen::VectorXi  IntVector;
typedef Eigen::VectorXd  DoubleVector;
typedef Eigen::VectorXcd ComplexVector;

typedef Eigen::MatrixXi  IntMatrix;
typedef Eigen::MatrixXd  DoubleMatrix;
typedef Eigen::MatrixXcd ComplexMatrix;

extern Getkw Input;
extern mrcpp::MultiResolutionAnalysis<3> *MRA;  //< Global MRA

} //namespace mrchem

