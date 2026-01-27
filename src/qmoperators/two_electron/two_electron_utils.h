#pragma once

#include "mrchem.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/qmfunction_fwd.h"
#include <unsupported/Eigen/CXX11/Tensor>

namespace mrchem {
/** @brief precomputes all the <ik|1/r12|jl> 2-electron integrals
 *
 *  @param[in] Phi input orbitals
 *
 * The ij/r12 potentials are first (pre)computed for all the orbitals j<i, then multiplied by kl
 */
Eigen::Tensor<std::complex<double>, 4> calc_2elintegrals(double prec, OrbitalVector &Phi);
/** @brief computes Int(phi_i^dag*phi_j/|r-r'|)
 *
 *  \param[in] phi_i orbital to be conjugated and multiplied by phi_j
 *  \param[in] phi_j orbital to be multiplied by phi_i^dag
 *  \param[in] rho is normally the density. It is used to screen the final output
 *  \param[out] V_ij is the resulting potential function
  *
 * Computes the product of complex conjugate of phi_i and phi_j,
 * then applies the Poisson operator. The result is given in V_ij.
 */
void ij_r12(double prec, Orbital rho, Orbital phi_i, Orbital phi_j, Orbital &V_ij) ;
}
