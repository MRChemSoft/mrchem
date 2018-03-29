#include <Eigen/Eigenvalues>

#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"

#include "GroundStateSolver.h"
#include "FockOperator.h"
#include "Orbital.h"

using namespace std;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

GroundStateSolver::GroundStateSolver(HelmholtzVector &h)
        : SCF(h),
          fMat_n(0),
          fOper_n(0),
          orbitals_n(0) {
}

GroundStateSolver::~GroundStateSolver() {
    if (this->fMat_n != 0) MSG_ERROR("Solver not properly cleared");
    if (this->fOper_n != 0) MSG_ERROR("Solver not properly cleared");
    if (this->orbitals_n != 0) MSG_ERROR("Solver not properly cleared");
}

/** Computes the Helmholtz argument for the all orbitals.
 *
 * Argument contains the potential operator acting on orbital i, and the sum
 * of all orbitals weighted by the Fock matrix. The effect of using inexact
 * Helmholtz operators are included in Lambda, which is a diagonal matrix
 * with the actual lambda parameters used in the Helmholtz operators
 * (input matrix M is assumed to be L-F).
 *
 * greenArg = \hat{V}orb_i + \sum_j (\Lambda_{ij}-F_{ij})orb_j
 */
OrbitalVector GroundStateSolver::setupHelmholtzArguments(FockOperator &fock,
                                                         const ComplexMatrix &M,
                                                         OrbitalVector &Phi,
                                                         bool adjoint,
                                                         bool clearFock) {
    Timer timer_tot;
    Printer::printHeader(0, "Setting up Helmholtz arguments");
    int oldprec = Printer::setPrecision(5);

    RankZeroTensorOperator &V = fock.potential();
    
    double coef = -1.0/(2.0*MATHCONST::pi);
    Timer timer_1;
    OrbitalVector part_1 = orbital::param_copy(Phi);
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) {
            part_1[i] = V(Phi[i]);
        }
    }
    timer_1.stop();
    Printer::printDouble(0, "Potential part", timer_1.getWallTime());

    if (clearFock) fock.clear();

    Timer timer_2;
    OrbitalVector part_2 = orbital::multiply(M, Phi);
    timer_2.stop();
    Printer::printDouble(0, "Matrix part", timer_2.getWallTime());

    OrbitalVector out = orbital::add(coef, part_1, coef, part_2);
    orbital::free(part_1);
    orbital::free(part_2);

    /*
    Timer timer_2;
    OrbitalVector orbVecChunk_i(0); //to store adresses of own i_orbs
    OrbitalVector rcvOrbs(0);       //to store adresses of received orbitals
    vector<int> orbsIx;             //to store own orbital indices
    int rcvOrbsIx[workOrbVecSize];  //to store received orbital indices

    //make vector with adresses of own orbitals
    int Ni = phi.size();
    int maxOrbPerMpi = Ni/mpiOrbSize + 1;//upper bound
    for (int ix = mpiOrbRank; ix < Ni; ix += mpiOrbSize) {
        orbVecChunk_i.push_back(phi.getOrbital(ix));//i orbitals
        orbsIx.push_back(ix);
    }

    for (int iter = 0; iter >= 0; iter++) {
        //get a new chunk from other processes
        orbVecChunk_i.getOrbVecChunk(orbsIx, rcvOrbs, rcvOrbsIx, Ni, iter);
        for (int i = mpiOrbRank; i < Ni; i += mpiOrbSize) {
            Orbital &phi_i = phi.getOrbital(i);

            vector<complex<double> > coefs;
            vector<Orbital *> orbs;

            int ix = i;
            for (int j = 0; j < rcvOrbs.size(); j++) {
                int jx = rcvOrbsIx[j];
                double coef = M(ix,jx);
                // Linear scaling screening inserted here
                // if (std::abs(coef) > MachineZero) {
                    Orbital &phi_j = rcvOrbs.getOrbital(j);
                    double norm_j = sqrt(phi_j.getSquareNorm());
                    if (norm_j > 0.01*getOrbitalPrecision()) {
                        coefs.push_back(coef);
                        orbs.push_back(&phi_j);
                    }
                //}
            }

            Orbital *tmp_i = new Orbital(phi_i);
            if (orbs.size() > 0) this->add(*tmp_i, coefs, orbs, false);

            this->add.inPlace(out[i], coef, *tmp_i);
            delete tmp_i;
        }
        rcvOrbs.clearVec(false);//reset to zero size orbital vector
    }
    orbVecChunk_i.clearVec(false);
    workOrbVec.clear();
    timer_2.stop();
    Printer::printDouble(0, "Matrix part", timer_2.getWallTime());
    */

    timer_tot.stop();
    Printer::printFooter(0, timer_tot, 2);
    Printer::setPrecision(oldprec);

    return out;
}

double GroundStateSolver::calcProperty() {
    Printer::printHeader(0, "Calculating SCF energy");
    Timer timer;

    ComplexMatrix &F = *this->fMat_n;
    FockOperator &fock = *this->fOper_n;
    OrbitalVector &Phi = *this->orbitals_n;

    SCFEnergy E = fock.trace(Phi, F);
    this->energy.push_back(E);

    timer.stop();
    int oldPrec = Printer::setPrecision(15);
    println(0, " Nuclear energy              " << setw(30) << E.getNuclearEnergy());
    println(0, " Electronic energy           " << setw(30) << E.getElectronicEnergy());
    Printer::printSeparator(0, '-');
    println(0, " Total energy                " << setw(30) << E.getTotalEnergy());
    Printer::printFooter(0, timer, 2);
    Printer::setPrecision(oldPrec);

    return E.getTotalEnergy();
}

double GroundStateSolver::calcPropertyError() const {
    int iter = this->property.size();
    return std::abs(getUpdate(this->property, iter, false));
}

void GroundStateSolver::printProperty() const {
    SCFEnergy scf_0, scf_1;
    int iter = this->energy.size();
    if (iter > 1) scf_0 = this->energy[iter - 2];
    if (iter > 0) scf_1 = this->energy[iter - 1];

    double T_0 = scf_0.getKineticEnergy();
    double T_1 = scf_1.getKineticEnergy();
    double V_0 = scf_0.getElectronNuclearEnergy();
    double V_1 = scf_1.getElectronNuclearEnergy();
    double J_0 = scf_0.getElectronElectronEnergy();
    double J_1 = scf_1.getElectronElectronEnergy();
    double K_0 = scf_0.getExchangeEnergy();
    double K_1 = scf_1.getExchangeEnergy();
    double XC_0 = scf_0.getExchangeCorrelationEnergy();
    double XC_1 = scf_1.getExchangeCorrelationEnergy();
    double E_0 = scf_0.getElectronicEnergy();
    double E_1 = scf_1.getElectronicEnergy();
    double N_0 = scf_0.getNuclearEnergy();
    double N_1 = scf_1.getNuclearEnergy();

    Printer::printHeader(0, "                    Energy                 Update      Done ");
    printUpdate(" Kinetic    ",  T_1,  T_1 -  T_0);
    printUpdate(" N-E        ",  V_1,  V_1 -  V_0);
    printUpdate(" Coulomb    ",  J_1,  J_1 -  J_0);
    printUpdate(" Exchange   ",  K_1,  K_1 -  K_0);
    printUpdate(" X-C        ", XC_1, XC_1 - XC_0);
    Printer::printSeparator(0, '-');
    printUpdate(" Electronic ",  E_1,  E_1 -  E_0);
    printUpdate(" Nuclear    ",  N_1,  N_1 -  N_0);
    Printer::printSeparator(0, '-');
    printUpdate(" Total      ",  E_1 + N_1, (E_1+N_1) - (E_0+N_0));
    Printer::printSeparator(0, '=');
}

/** Minimize the spatial extension of orbitals, by a transformation of orbitals
 *
 * Minimizes \f$  \sum_{i=1,N}\langle i| {\bf R^2}  | i \rangle - \langle i| {\bf R}| i \rangle^2 \f$
 *	which is equivalent to maximizing \f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 *
 * The resulting transformation includes the orthonormalization of the orbitals.
 * For details see the tex documentation in doc directory
 */
void GroundStateSolver::localize(double prec, FockOperator &fock, ComplexMatrix &F, OrbitalVector &Phi) {
    NOT_IMPLEMENTED_ABORT;
    /*
    Printer::printHeader(0, "Localizing orbitals");
    Timer timer;

    MatrixXd U;
    int n_it = 0;
    if (Phi.size() > 1) {
        Timer rr_t;
        RR rr(this->orbPrec[0], Phi);
        n_it = rr.maximize();//compute total U, rotation matrix
        rr_t.stop();
        Printer::printDouble(0, "Computing Foster-Boys matrix", rr_t.getWallTime());
        if (n_it > 0) {
            println(0, " Converged after iteration   " << setw(30) << n_it);
            U = rr.getTotalU().transpose();
        } else {
            println(0, " Foster-Boys localization did not converge!");
        }
    } else {
        println(0, " Cannot localize less than two orbitals");
    }

    if (n_it <= 0) {
        Timer orth_t;
        U = calcOrthonormalizationMatrix(Phi);
        orth_t.stop();
        Printer::printDouble(0, "Computing Lowdin matrix", orth_t.getWallTime());
    }

    Timer rot_t;
    F = U*F*U.transpose();
    fock.rotate(U);
    this->add.rotate(Phi, U);
    rot_t.stop();
    Printer::printDouble(0, "Rotating orbitals", rot_t.getWallTime());

    timer.stop();
    Printer::printFooter(0, timer, 2);
    */
}

/** Perform the orbital rotation that diagonalizes the Fock matrix
 *
 * This operation includes the orthonormalization using the overlap matrix.
 */
void GroundStateSolver::diagonalize(double prec, FockOperator &fock, ComplexMatrix &F, OrbitalVector &Phi) {
    Printer::printHeader(0, "Digonalizing Fock matrix");
    Timer timer;

    Timer orth_t;
    ComplexMatrix S_m12 = calcOrthonormalizationMatrix(Phi);
    F = S_m12.transpose()*F*S_m12;
    orth_t.stop();
    Printer::printDouble(0, "Computing Lowdin matrix", orth_t.getWallTime());

    Timer diag_t;
    ComplexMatrix U = ComplexMatrix::Zero(F.rows(), F.cols());
    int np = orbital::size_paired(Phi);
    int na = orbital::size_alpha(Phi);
    int nb = orbital::size_beta(Phi);
    if (np > 0) diagonalizeBlock(F, U, 0,       np);
    if (na > 0) diagonalizeBlock(F, U, np,      na);
    if (nb > 0) diagonalizeBlock(F, U, np + na, nb);
    U = U * S_m12;
    diag_t.stop();
    Printer::printDouble(0, "Diagonalizing matrix", diag_t.getWallTime());

    Timer rot_t;
    fock.rotate(U);
    OrbitalVector Psi = orbital::multiply(U, Phi, prec);
    orbital::free(Phi);
    Phi = Psi;
    rot_t.stop();
    Printer::printDouble(0, "Rotating orbitals", rot_t.getWallTime());

    timer.stop();
    Printer::printFooter(0, timer, 2);
}


void GroundStateSolver::diagonalizeBlock(ComplexMatrix &M, ComplexMatrix &U, int nstart, int nsize) {
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(nsize);
    es.compute(M.block(nstart, nstart, nsize, nsize));
    ComplexMatrix tmp = es.eigenvectors();
    U.block(nstart, nstart, nsize, nsize) = tmp.transpose();
    M.block(nstart, nstart, nsize, nsize) = es.eigenvalues().asDiagonal();
}

void GroundStateSolver::orthonormalize(double prec, FockOperator &fock, ComplexMatrix &F, OrbitalVector &Phi) {
    ComplexMatrix U = calcOrthonormalizationMatrix(Phi);
    F = U*F*U.transpose();
    fock.rotate(U);
    OrbitalVector Psi = orbital::multiply(U, Phi, prec);
    orbital::free(Phi);
    Phi = Psi;
}

ComplexMatrix GroundStateSolver::calcOrthonormalizationMatrix(OrbitalVector &Phi) {
    Timer timer;
    printout(1, "Calculating orthonormalization matrix            ");

    ComplexMatrix S_tilde = orbital::calc_overlap_matrix(Phi);

    Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(S_tilde.cols());
    es.compute(S_tilde);

    DoubleMatrix A = es.eigenvalues().asDiagonal();
    for (int i = 0; i < A.cols(); i++) {
        A(i,i) = std::pow(A(i,i), -1.0/2.0);
    }
    ComplexMatrix B = es.eigenvectors();
    ComplexMatrix U = B*A*B.transpose();

    timer.stop();
    println(1, timer.getWallTime());
    return U;
}

/** Compute the position matrix <i|R_x|j>,<i|R_y|j>,<i|R_z|j>
 */
RR::RR(double prec, OrbitalVector &Phi) {
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
double RR::functional() {
    //    //s1 is what should be maximized (i.e. the sum of <i R i>^2)
    double s1=0.0;
    for (int dim=0; dim<3; dim++) {
        for (int j=0; j<N; j++) {
            s1+=r_i(j,j+dim*N)*r_i(j,j+dim*N);
        }
    }
    return s1;
}

/** Make gradient vector for the RR case
 */
double RR::make_gradient() {
    double norm = 0.0;
    for (int i=0; i<N2h; i++) gradient(i)=0.0 ;
    for (int dim=0; dim<3; dim++) {
        int ij=0;
        for (int j=0; j<N; j++) {
            for (int i=0; i<j; i++) {
                gradient(ij)+=4.0*r_i(i,j+dim*N)*(r_i(i,i+dim*N)-r_i(j,j+dim*N));
                ij++;
            }
        }
    }
    for (int ij=0; ij<N2h; ij++) {
        norm += gradient(ij)*gradient(ij);
    }
    return std::sqrt(norm);
}


/** Make Hessian matrix for the RR case
 */
double RR::make_hessian() {
    double djk,djl,dik,dil;
    for (int j=0; j<N2h; j++){
        for (int i=0; i<N2h; i++){
            hessian(i,j)=0.0 ;
        }
    }
    for (int dim=0; dim<3; dim++) {
        int kl=0;
        for (int l=0; l<N; l++) {
            for (int k=0; k<l; k++) {
                int ij=0;
                for (int j=0; j<N; j++) {
                    for (int i=0; i<j; i++) {
                        djk = j==k ? 1.0: 0.0;
                        djl = j==l ? 1.0: 0.0;
                        dik = i==k ? 1.0: 0.0;
                        dil = i==l ? 1.0: 0.0;

                        hessian(ij,kl)+=2.0*(
                                    djk*r_i(i,i+dim*N)*r_i(l,i+dim*N)
                                    -djl*r_i(i,i+dim*N)*r_i(k,i+dim*N)
                                    -dik*r_i(j,j+dim*N)*r_i(l,j+dim*N)
                                    +dil*r_i(j,j+dim*N)*r_i(k,j+dim*N)
                                    -2*dil*r_i(i,i+dim*N)*r_i(k,j+dim*N)
                                    +2.0*dik*r_i(i,i+dim*N)*r_i(l,j+dim*N)
                                    +2*djl*r_i(j,j+dim*N)*r_i(k,i+dim*N)
                                    -2.0*djk*r_i(j,j+dim*N)*r_i(l,i+dim*N)
                                    +djk*r_i(l,l+dim*N)*r_i(i,l+dim*N)
                                    -dik*r_i(l,l+dim*N)*r_i(j,l+dim*N)
                                    -djl*r_i(k,k+dim*N)*r_i(i,k+dim*N)
                                    +dil*r_i(k,k+dim*N)*r_i(j,k+dim*N)
                                    -4*(dil-dik-djl+djk)*r_i(i,j+dim*N)*r_i(k,l+dim*N));

                        ij++;
                    }
                }
                kl++;
            }
        }
    }

    return 0;
}

/** Given the step matrix, update the rotation matrix and the R matrix
 */
void RR::do_step(DoubleVector step) {
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

