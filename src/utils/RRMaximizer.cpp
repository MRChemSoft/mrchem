#include "RRMaximizer.h"

namespace mrchem {

/** Compute the position matrix <i|R_x|j>,<i|R_y|j>,<i|R_z|j>
 */
RRMaximizer::RRMaximizer(double prec, OrbitalVector &Phi) {
    NOT_IMPLEMENTED_ABORT;
    /*
    N = Phi.size();
    if (N < 2) MSG_ERROR("Cannot localize less than two orbitals");
    total_U = MatrixXd::Identity(N,N);
    N2h = N*(N-1)/2;
    gradient = VectorXd(N2h);
    hessian = MatrixXd(N2h, N2h);
    r_i_orig = MatrixXd::Zero(N,3*N);
    r_i = MatrixXd(N,3*N);

    //Make R matrix
    PositionOperator r;
    r.setup(prec);

    RankZeroTensorOperator &r_x = r[0];
    RankZeroTensorOperator &r_y = r[1];
    RankZeroTensorOperator &r_z = r[2];

#ifdef HAVE_MPI

    OrbitalVector orbVecChunk_i(0); //to store adresses of own i_orbs
    OrbitalVector rcvOrbs(0);       //to store adresses of received orbitals
    vector<int> orbsIx;             //to store own orbital indices
    int rcvOrbsIx[workOrbVecSize];  //to store received orbital indices

    //make vector with adresses of own orbitals
    for (int Ix = mpiOrbRank; Ix < N; Ix += mpiOrbSize) {
        orbVecChunk_i.push_back(Phi.getOrbital(Ix));//i orbitals
        orbsIx.push_back(Ix);
    }

    for (int iter = 0; iter >= 0; iter++) {
        //get a new chunk from other processes
        orbVecChunk_i.getOrbVecChunk_sym(orbsIx, rcvOrbs, rcvOrbsIx, N, iter);
        for (int i = 0; i<orbVecChunk_i.size(); i++){
            Orbital &phi_i = orbVecChunk_i.getOrbital(i);
            int spin_i = phi_i.getSpin();
            for (int j = 0; j < rcvOrbs.size(); j++) {
                Orbital &phi_j = rcvOrbs.getOrbital(j);
                int spin_j = phi_j.getSpin();
                if (spin_i != spin_j) {
                    MSG_ERROR("Spins must be separated before localization");
                }
                //NOTE: the "if" should not be necessary, but since outside the required precision
                //r_x(phi_i, phi_j) != r_x(phi_j, phi_i), we prefer to have consistent results for
                //different mpiOrbSize
                if(rcvOrbsIx[j]<=orbsIx[i]){
                    r_i_orig(orbsIx[i],rcvOrbsIx[j]) = r_x(phi_i, phi_j);
                    r_i_orig(orbsIx[i],rcvOrbsIx[j]+N) = r_y(phi_i, phi_j);
                    r_i_orig(orbsIx[i],rcvOrbsIx[j]+2*N) = r_z(phi_i, phi_j);
                    r_i_orig(rcvOrbsIx[j],orbsIx[i]) = r_i_orig(orbsIx[i],rcvOrbsIx[j]);
                    r_i_orig(rcvOrbsIx[j],orbsIx[i]+N) =  r_i_orig(orbsIx[i],rcvOrbsIx[j]+N);
                    r_i_orig(rcvOrbsIx[j],orbsIx[i]+2*N) = r_i_orig(orbsIx[i],rcvOrbsIx[j]+2*N);
                }else{
                    if(rcvOrbsIx[j]%mpiOrbSize != mpiOrbRank){ //need only compute j<=i in own block
                        r_i_orig(rcvOrbsIx[j],orbsIx[i]) = r_x(phi_j, phi_i);
                        r_i_orig(rcvOrbsIx[j],orbsIx[i]+N) =  r_y(phi_j, phi_i);
                        r_i_orig(rcvOrbsIx[j],orbsIx[i]+2*N) = r_z(phi_j, phi_i);
                        r_i_orig(orbsIx[i],rcvOrbsIx[j]) = r_i_orig(rcvOrbsIx[j],orbsIx[i]);
                        r_i_orig(orbsIx[i],rcvOrbsIx[j]+N) = r_i_orig(rcvOrbsIx[j],orbsIx[i]+N) ;
                        r_i_orig(orbsIx[i],rcvOrbsIx[j]+2*N) = r_i_orig(rcvOrbsIx[j],orbsIx[i]+2*N);
                    }
                }
            }
        }
        rcvOrbs.clearVec(false);
    }
    orbVecChunk_i.clearVec(false);
    workOrbVec.clear();
    //combine results from all processes
    MPI_Allreduce(MPI_IN_PLACE, &r_i_orig(0,0), N*3*N,
                  MPI_DOUBLE, MPI_SUM, mpiCommOrb);

#else

    for (int i = 0; i < N; i++) {
        Orbital &phi_i = phi.getOrbital(i);
        int spin_i = phi_i.getSpin();
        for (int j = 0; j <= i; j++) {
            Orbital &phi_j =  phi.getOrbital(j);
            int spin_j = phi_j.getSpin();
            if (spin_i != spin_j) {
                MSG_ERROR("Spins must be separated before localization");
            }
            r_i_orig(i,j) = r_x(phi_i, phi_j);
            r_i_orig(i,j+N) = r_y(phi_i, phi_j);
            r_i_orig(i,j+2*N) = r_z(phi_i, phi_j);
            r_i_orig(j,i) = r_i_orig(i,j);
            r_i_orig(j,i+N) = r_i_orig(i,j+N);
            r_i_orig(j,i+2*N) = r_i_orig(i,j+2*N);
        }
    }
#endif
    r_x.clear();
    r_y.clear();
    r_z.clear();

    //rotate R matrices into orthonormal basis
    IdentityOperator I;
    I.setup(prec);
    MatrixXd S_tilde = I(phi, phi);
    MatrixXd S_tilde_m12 = MathUtils::hermitianMatrixPow(S_tilde, -1.0/2.0);
    I.clear();

    total_U=S_tilde_m12*total_U;
    MatrixXd R(N, N);
    for(int dim=0; dim<3; dim++){
        for (int j=0; j<N; j++) {
            for (int i=0; i<=j; i++) {
                R(i,j)=r_i_orig(i,j+dim*N);
                R(j,i)=r_i_orig(i,j+dim*N);//Enforce symmetry
            }
        }
        R=total_U.transpose()*R*total_U;
        for (int j=0; j<N; j++) {
            for (int i=0; i<=j; i++) {
                r_i(i,j+dim*N)=R(i,j);
                r_i(j,i+dim*N)=R(i,j);//Enforce symmetry
            }
        }
    }
    */
}

/** compute the value of
 * f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 */
double RRMaximizer::functional() const {
    //s1 is what should be maximized (i.e. the sum of <i R i>^2)
    double s1 = 0.0;
    for (int d = 0; d < 3; d++) {
        for (int j = 0; j < this->N; j++) {
            s1 += this->r_i(j,j+d*this->N)*this->r_i(j,j+d*this->N);
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
                this->gradient(ij) += 4.0*r_i(i,j+d*N)*(r_i(i,i+d*N)-r_i(j,j+d*N));
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
void RRMaximizer::do_step(DoubleVector &step) {
    NOT_IMPLEMENTED_ABORT;
    /*
    MatrixXd A(N,N);
    //define rotation U=exp(-A), A real antisymmetric, from step
    int ij=0;
    for (int j=0; j<N; j++) {
        for (int i=0; i<j; i++) {
            A(i,j)=step(ij);
            A(j,i)=-A(i,j);
            ij++;
        }
        A(j,j)=0.0;
    }

    //calculate U=exp(-A) by diagonalization and U=Vexp(id)Vt with VdVt=iA
    //could also sum the term in the expansion if A is small
    total_U*=MathUtils::SkewMatrixExp(A);

    //rotate the original r matrix with total U
    MatrixXd r(N, N);
    for(int dim=0; dim<3; dim++){
        for (int j=0; j<N; j++) {
            for (int i=0; i<N; i++) {
                r(i,j)=r_i_orig(i,j+dim*N);
            }
        }
        r=total_U.transpose()*r*total_U;
        for (int j=0; j<N; j++) {
            for (int i=0; i<N; i++) {
                r_i(i,j+dim*N)=r(i,j);
            }
        }
    }
    */
}

} //namespace mrchem

