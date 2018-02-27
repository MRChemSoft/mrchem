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

#include "MRCPP/MWFunctions"

namespace mrchem {

typedef std::complex<double> ComplexDouble;

extern mrcpp::MultiResolutionAnalysis<3> *MRA;  //< Default MRA

} //namespace mrchem

