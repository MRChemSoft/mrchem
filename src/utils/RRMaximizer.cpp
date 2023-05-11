/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "MRCPP/Printer"
#include <unsupported/Eigen/MatrixFunctions> // faster exponential of matrices

#include "utils/RRMaximizer.h"
#include "utils/math_utils.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/OrbitalIterator.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/PositionOperator.h"

namespace mrchem {

/** Compute the position matrix <i|R_x|j>,<i|R_y|j>,<i|R_z|j>
 */
RRMaximizer::RRMaximizer(double prec, OrbitalVector &Phi) {
    this->N = Phi.size();
    if (this->N < 2) MSG_ERROR("Cannot localize less than two orbitals");

    this->total_U = DoubleMatrix::Identity(this->N, this->N);
    this->N2h = this->N * (this->N - 1) / 2;
    this->gradient = DoubleVector(this->N2h);
    this->r_i_orig = DoubleMatrix::Zero(this->N, 3 * this->N);
    this->r_i = DoubleMatrix(this->N, 3 * this->N);

    // Make R matrix
    PositionOperator r;
    r.setup(prec);

    RankZeroOperator &r_x = r[0];
    RankZeroOperator &r_y = r[1];
    RankZeroOperator &r_z = r[2];

    ComplexMatrix R_x = ComplexMatrix::Zero(Phi.size(), Phi.size());
    ComplexMatrix R_y = ComplexMatrix::Zero(Phi.size(), Phi.size());
    ComplexMatrix R_z = ComplexMatrix::Zero(Phi.size(), Phi.size());

    OrbitalVector xPhi_Vec = r_x(Phi);
    OrbitalVector yPhi_Vec = r_y(Phi);
    OrbitalVector zPhi_Vec = r_z(Phi);

    R_x = orbital::calc_overlap_matrix(Phi, xPhi_Vec);
    R_y = orbital::calc_overlap_matrix(Phi, yPhi_Vec);
    R_z = orbital::calc_overlap_matrix(Phi, zPhi_Vec);

    for (int i = 0; i < this->N; i++) {
        for (int j = 0; j <= i; j++) {
            if (Phi[i].spin() != Phi[j].spin()) MSG_ERROR("Spins must be separated before localization");
            this->r_i_orig(i, j) = R_x(i, j).real();
            this->r_i_orig(i, j + this->N) = R_y(i, j).real();
            this->r_i_orig(i, j + 2 * this->N) = R_z(i, j).real();

            this->r_i_orig(j, i) = this->r_i_orig(i, j);
            this->r_i_orig(j, i + this->N) = this->r_i_orig(i, j + this->N);
            this->r_i_orig(j, i + 2 * this->N) = this->r_i_orig(i, j + 2 * this->N);
        }
    }
    r_x.clear();
    r_y.clear();
    r_z.clear();

    // rotate R matrices into orthonormal basis
    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);

    this->total_U = S_m12.real() * this->total_U;

    DoubleMatrix R(this->N, this->N);
    for (int d = 0; d < 3; d++) {
        for (int j = 0; j < this->N; j++) {
            for (int i = 0; i <= j; i++) {
                R(i, j) = this->r_i_orig(i, j + d * this->N);
                R(j, i) = this->r_i_orig(i, j + d * this->N); // Enforce symmetry
            }
        }
        R = this->total_U.transpose() * R * this->total_U;
        for (int j = 0; j < this->N; j++) {
            for (int i = 0; i <= j; i++) {
                this->r_i(i, j + d * this->N) = R(i, j);
                this->r_i(j, i + d * this->N) = R(i, j); // Enforce symmetry
            }
        }
    }
}

/** compute the value of
 * f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 */
double RRMaximizer::functional() const {
    // s1 is what should be maximized (i.e. the sum of <i R i>^2)
    double s1 = 0.0;
    for (int d = 0; d < 3; d++) {
        for (int j = 0; j < this->N; j++) {
            double r_jj = this->r_i(j, j + d * this->N);
            s1 += r_jj * r_jj;
        }
    }
    return s1;
}

/** Make gradient vector for the RRMaximizer case
 */
double RRMaximizer::make_gradient() {
    this->gradient = DoubleVector::Zero(this->N2h);

    for (int d = 0; d < 3; d++) {
        int ij = 0;
        for (int j = 0; j < this->N; j++) {
            for (int i = 0; i < j; i++) {
                double r_ij = this->r_i(i, j + d * this->N);
                double r_ii = this->r_i(i, i + d * this->N);
                double r_jj = this->r_i(j, j + d * this->N);
                this->gradient(ij) += 4.0 * r_ij * (r_ii - r_jj);
                ij++;
            }
        }
    }

    double norm = 0.0;
    for (int ij = 0; ij < this->N2h; ij++) { norm += this->gradient(ij) * this->gradient(ij); }
    return std::sqrt(norm);
}

/** Make Hessian matrix for the RRMaximizer case
 */
// clang-format off
double RRMaximizer::make_hessian() {
    MSG_ABORT("Hessian is not explicitely constructed");

    hessian = DoubleMatrix::Zero(this->N2h, this->N2h);

    double djk, djl, dik, dil;
    for (int d = 0; d < 3; d++) {
        int kl = 0;
        for (int l = 0; l < this->N; l++) {
            for (int k = 0; k < l; k++) {
                int ij = 0;
                for (int j = 0; j < this->N; j++) {
                    for (int i = 0; i < j; i++) {
                        djk = (j == k) ? 1.0 : 0.0;
                        djl = (j == l) ? 1.0 : 0.0;
                        dik = (i == k) ? 1.0 : 0.0;
                        dil = (i == l) ? 1.0 : 0.0;

                        hessian(ij,kl) += 2.0*(
                                         djk*r_i(i,i+d*N)*r_i(l,i+d*N)
                                        -djl*r_i(i,i+d*N)*r_i(k,i+d*N)
                                        -dik*r_i(j,j+d*N)*r_i(l,j+d*N)
                                        +dil*r_i(j,j+d*N)*r_i(k,j+d*N)
                                    -2.0*dil*r_i(i,i+d*N)*r_i(k,j+d*N)
                                    +2.0*dik*r_i(i,i+d*N)*r_i(l,j+d*N)
                                    +2.0*djl*r_i(j,j+d*N)*r_i(k,i+d*N)
                                    -2.0*djk*r_i(j,j+d*N)*r_i(l,i+d*N)
                                        +djk*r_i(l,l+d*N)*r_i(i,l+d*N)
                                        -dik*r_i(l,l+d*N)*r_i(j,l+d*N)
                                        -djl*r_i(k,k+d*N)*r_i(i,k+d*N)
                                        +dil*r_i(k,k+d*N)*r_i(j,k+d*N)
                      -4.0*(dil-dik-djl+djk)*r_i(i,j+d*N)*r_i(k,l+d*N));

                        ij++;
                    }
                }
                kl++;
            }
        }
    }
    return 0.0;
}
// clang-format on

/** Multiply a vector with Hessian for the RRMaximizer case.
 *  Does not store the Hessian values.
 */
// clang-format off
DoubleMatrix unflattenSkew(const DoubleVector &v, int N) {
  DoubleMatrix M(N,N);
  int k = 0;
  for (int i = 0; i < N; i++) {
    M(i,i) = 0;
    for (int j = 0; j < i; j++) {
      M(i,j) = v(k++);
      M(j,i) = -M(i,j);
    }
  }
  return M;
}

DoubleVector flattenSkew(const DoubleMatrix &M, int N) {
  DoubleVector v( N*(N-1)/2 );
  int k = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < i; j++) {
      v(k++) = M(i,j);
    }
  }
  return v;
}

void RRMaximizer::multiply_hessian(DoubleVector &vec, DoubleVector &Hv) {

  DoubleMatrix S = unflattenSkew(vec, this->N);
  DoubleMatrix HS(this->N, this->N);
  HS.setZero();

  // fast and simple method by Johan S. Wind:
  for (int d = 0; d < 3; d++) {
    // D = diag(R);
    // RS = R@S;
    // return (4*D@S - 2*S@D - 8*diag(RS))@R - 2*D@RS;

    DoubleMatrix R = r_i.block(0,N*d,this->N,this->N);
    DoubleMatrix RS = R * S;

    auto&D = R.diagonal().asDiagonal();
    auto&RS_diag = RS.diagonal().asDiagonal();

    HS += (4*D*S - 2*S*D) * R - 2*D*RS - 8*RS_diag*R;
  }
  Hv = flattenSkew(HS-HS.transpose(), this->N);

}
// clang-format on

/** Returns Hessian value for the RRMaximizer case
 */
// clang-format off
double RRMaximizer::get_hessian(int ij, int kl) {
    double hessian_val = 0.0;

    double jr = 0.5*(std::sqrt(9.0+8.0*ij)+1.0);
    auto j = (int) (jr-1.0E-9);//the last jr is exactly the integer *over* the wanted j, if something is not subtracted
    int i = ij-(j*(j-1))/2;

    if(i>=j or i<0 or j<0)std::cout<<" ij ERROR "<<i<<" "<<j<<" "<<ij<<std::endl;
    //std::cout<<"i "<<i<<" j"<<j<<" ij"<<ij<<std::endl;
    double lr = 0.5*(std::sqrt(9.0+8.0*kl)+1.0);
    auto l = (int) (lr-1.0E-9);
    int k = kl-(l*(l-1))/2;

    if(k>=l or k<0 or l<0)std::cout<<" kl ERROR "<<k<<" "<<l<<" "<<kl<<std::endl;

    double djk, djl, dik, dil;
    for (int d = 0; d < 3; d++) {
        djk = (j == k) ? 1.0 : 0.0;
        djl = (j == l) ? 1.0 : 0.0;
        dik = (i == k) ? 1.0 : 0.0;
        dil = (i == l) ? 1.0 : 0.0;

        hessian_val += 2.0*(
                    djk*r_i(i,i+d*N)*r_i(l,i+d*N)
                    -djl*r_i(i,i+d*N)*r_i(k,i+d*N)
                    -dik*r_i(j,j+d*N)*r_i(l,j+d*N)
                    +dil*r_i(j,j+d*N)*r_i(k,j+d*N)
                    -2.0*dil*r_i(i,i+d*N)*r_i(k,j+d*N)
                    +2.0*dik*r_i(i,i+d*N)*r_i(l,j+d*N)
                    +2.0*djl*r_i(j,j+d*N)*r_i(k,i+d*N)
                    -2.0*djk*r_i(j,j+d*N)*r_i(l,i+d*N)
                    +djk*r_i(l,l+d*N)*r_i(i,l+d*N)
                    -dik*r_i(l,l+d*N)*r_i(j,l+d*N)
                    -djl*r_i(k,k+d*N)*r_i(i,k+d*N)
                    +dil*r_i(k,k+d*N)*r_i(j,k+d*N)
                    -4.0*(dil-dik-djl+djk)*r_i(i,j+d*N)*r_i(k,l+d*N));
    }
    return hessian_val;
}
// clang-format on

/** Given the step matrix, update the rotation matrix and the R matrix
 */
void RRMaximizer::do_step(const DoubleVector &step) {
    DoubleMatrix A(this->N, this->N);
    // define rotation U=exp(-A), A real antisymmetric, from step
    int ij = 0;
    for (int j = 0; j < this->N; j++) {
        for (int i = 0; i < j; i++) {
            A(i, j) = step(ij);
            A(j, i) = -step(ij);
            ij++;
        }
        A(j, j) = 0.0;
    }

    // calculate U=exp(-A) by diagonalization and U=Vexp(id)Vt with VdVt=iA
    // could also sum the terms in the expansion
    this->total_U *= A.exp().transpose(); // transpose because we want exp(-A)

    // rotate the original r matrix with total U
    DoubleMatrix r(this->N, this->N);
    for (int d = 0; d < 3; d++) {
        for (int j = 0; j < this->N; j++) {
            for (int i = 0; i < this->N; i++) { r(i, j) = this->r_i_orig(i, j + d * this->N); }
        }
        r = this->total_U.transpose() * r * this->total_U;
        for (int j = 0; j < this->N; j++) {
            for (int i = 0; i < this->N; i++) { this->r_i(i, j + d * this->N) = r(i, j); }
        }
    }
}

} // namespace mrchem
