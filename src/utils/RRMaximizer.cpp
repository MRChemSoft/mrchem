#include "MRCPP/Printer"

#include "utils/RRMaximizer.h"
#include "utils/math_utils.h"

#include "PositionOperator.h"

namespace mrchem {

/** Compute the position matrix <i|R_x|j>,<i|R_y|j>,<i|R_z|j>
 */
RRMaximizer::RRMaximizer(double prec, OrbitalVector &Phi) {
    this->N = Phi.size();
    if (this->N < 2) MSG_ERROR("Cannot localize less than two orbitals");

    this->total_U = DoubleMatrix::Identity(this->N,this->N);
    this->N2h = this->N*(this->N-1)/2;
    this->gradient = DoubleVector(this->N2h);
    this->hessian = DoubleMatrix(this->N2h, this->N2h);
    this->r_i_orig = DoubleMatrix::Zero(this->N,3*this->N);
    this->r_i = DoubleMatrix(this->N,3*this->N);

    //Make R matrix
    PositionOperator r;
    r.setup(prec);

    RankZeroTensorOperator &r_x = r[0];
    RankZeroTensorOperator &r_y = r[1];
    RankZeroTensorOperator &r_z = r[2];

    ComplexMatrix R_x = r_x(Phi, Phi);
    ComplexMatrix R_y = r_y(Phi, Phi);
    ComplexMatrix R_z = r_z(Phi, Phi);

    for (int i = 0; i < this->N; i++) {
        for (int j = 0; j <= i; j++) {
            if (Phi[i].spin() != Phi[j].spin()) MSG_ERROR("Spins must be separated before localization");
            this->r_i_orig(i,j          ) = R_x(i,j).real();
            this->r_i_orig(i,j+  this->N) = R_y(i,j).real();
            this->r_i_orig(i,j+2*this->N) = R_z(i,j).real();

            this->r_i_orig(j,i          ) = this->r_i_orig(i,j          );
            this->r_i_orig(j,i+  this->N) = this->r_i_orig(i,j+  this->N);
            this->r_i_orig(j,i+2*this->N) = this->r_i_orig(i,j+2*this->N);
        }
    }
    r_x.clear();
    r_y.clear();
    r_z.clear();

    //rotate R matrices into orthonormal basis
    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);

    this->total_U = S_m12.real()*this->total_U;
    DoubleMatrix R(this->N, this->N);
    for(int d = 0; d < 3; d++){
        for (int j = 0; j < this->N; j++) {
            for (int i = 0; i <= j; i++) {
                R(i,j) = this->r_i_orig(i,j+d*this->N);
                R(j,i) = this->r_i_orig(i,j+d*this->N);//Enforce symmetry
            }
        }
        R = this->total_U.transpose()*R*this->total_U;
        for (int j = 0; j < this->N; j++) {
            for (int i = 0; i <= j; i++) {
                this->r_i(i,j+d*this->N)=R(i,j);
                this->r_i(j,i+d*this->N)=R(i,j);//Enforce symmetry
            }
        }
    }
}

/** compute the value of
 * f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 */
double RRMaximizer::functional() const {
    //s1 is what should be maximized (i.e. the sum of <i R i>^2)
    double s1 = 0.0;
    for (int d = 0; d < 3; d++) {
        for (int j = 0; j < this->N; j++) {
            double r_jj = this->r_i(j,j+d*this->N);
            s1 += r_jj*r_jj;
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
                double r_ij = this->r_i(i,j+d*this->N);
                double r_ii = this->r_i(i,i+d*this->N);
                double r_jj = this->r_i(j,j+d*this->N);
                this->gradient(ij) += 4.0*r_ij*(r_ii - r_jj);
                ij++;
            }
        }
    }

    double norm = 0.0;
    for (int ij = 0; ij < this->N2h; ij++) {
        norm += this->gradient(ij)*this->gradient(ij);
    }
    return std::sqrt(norm);
}

/** Make Hessian matrix for the RRMaximizer case
 */
double RRMaximizer::make_hessian() {
    this->hessian = DoubleMatrix::Zero(this->N2h, this->N2h);

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

                        this->hessian(ij,kl) += 2.0*(
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

/** Given the step matrix, update the rotation matrix and the R matrix
 */
void RRMaximizer::do_step(const DoubleVector &step) {
    DoubleMatrix A(this->N, this->N);
    //define rotation U=exp(-A), A real antisymmetric, from step
    int ij = 0;
    for (int j = 0; j < this->N; j++) {
        for (int i = 0; i < j; i++) {
            A(i,j) =  step(ij);
            A(j,i) = -step(ij);
            ij++;
        }
        A(j,j) = 0.0;
    }

    //calculate U=exp(-A) by diagonalization and U=Vexp(id)Vt with VdVt=iA
    //could also sum the term in the expansion if A is small
    this->total_U *= math_utils::skew_matrix_exp(A);

    //rotate the original r matrix with total U
    DoubleMatrix r(this->N, this->N);
    for (int d = 0; d < 3; d++){
        for (int j = 0; j < this->N; j++) {
            for (int i = 0; i < this->N; i++) {
                r(i,j) = this->r_i_orig(i,j+d*this->N);
            }
        }
        r = this->total_U.transpose()*r*this->total_U;
        for (int j = 0; j < this->N; j++) {
            for (int i = 0; i < this->N; i++) {
                this->r_i(i,j+d*this->N) = r(i,j);
            }
        }
    }
}

} //namespace mrchem

