#pragma once

#include "mrchem.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/qmfunction_fwd.h"
#include <unsupported/Eigen/CXX11/Tensor>

namespace mrchem {
Eigen::Tensor<std::complex<double>, 4> calc_2elintegrals(double prec, OrbitalVector &Phi);
void ij_r12(double prec, Orbital rho, Orbital phi_i, Orbital phi_j, Orbital &V_ij) ;
}
