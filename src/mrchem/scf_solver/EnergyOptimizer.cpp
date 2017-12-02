#include "EnergyOptimizer.h"
#include "OrbitalVector.h"
#include "FockOperator.h"
#include "NuclearPotential.h"
#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "XCOperator.h"
#include "HelmholtzOperatorSet.h"
#include "IdentityOperator.h"
#include "eigen_disable_warnings.h"

using namespace std;
using namespace Eigen;

extern OrbitalVector workOrbVec;

EnergyOptimizer::EnergyOptimizer(HelmholtzOperatorSet &h)
    : GroundStateSolver(h),
      fOper_np1(0) {
}

EnergyOptimizer::~EnergyOptimizer() {
    if (this->fOper_np1 != 0) MSG_ERROR("Solver not properly cleared");
}

void EnergyOptimizer::setup(FockOperator &fock, OrbitalVector &phi, MatrixXd &F,
                            FockOperator &fock_np1, OrbitalVector &phi_np1) {
    this->fMat_n = &F;
    this->fOper_n = &fock;
    this->fOper_np1 = &fock_np1;

    this->orbitals_n = &phi;
    this->orbitals_np1 = &phi_np1;
    this->dOrbitals_n = new OrbitalVector(phi);
}

void EnergyOptimizer::clear() {
    this->fMat_n = 0;
    this->fOper_n = 0;
    this->fOper_np1 = 0;

    delete this->dOrbitals_n;

    this->orbitals_n = 0;
    this->orbitals_np1 = 0;
    this->dOrbitals_n = 0;

    resetPrecision();
}

bool EnergyOptimizer::optimize() {
    MatrixXd &F_n = *this->fMat_n;
    FockOperator &fock = *this->fOper_n;
    OrbitalVector &phi_n = *this->orbitals_n;
    OrbitalVector &phi_np1 = *this->orbitals_np1;
    OrbitalVector &dPhi_n = *this->dOrbitals_n;
    HelmholtzOperatorSet &H = *this->helmholtz;

    double err_o = phi_n.getErrors().maxCoeff();
    double err_t = 1.0;
    double err_p = 1.0;

    this->add.setPrecision(this->orbPrec[2]/10.0);
    if (this->canonical) {
        diagonalize(fock, F_n, phi_n);
    } else {
        localize(fock, F_n, phi_n);
    }

    int nIter = 0;
    bool converged = false;
    while(nIter++ < this->maxIter or this->maxIter < 0) {
        // Initialize SCF cycle
        Timer timer;
        printCycle(nIter);
        adjustPrecision(err_o);

        double orb_prec = getOrbitalPrecision();

        // Compute electronic energy
        fock.setup(orb_prec);
        double E = calcProperty();
        this->property.push_back(E);

        // Setup Helmholtz operators and argument
        H.setup(orb_prec, F_n.diagonal());
        MatrixXd L_n = H.getLambda().asDiagonal();
        OrbitalVector *args_n = setupHelmholtzArguments(fock, L_n-F_n, phi_n);

        // Apply Helmholtz operators
        H(phi_np1, *args_n);
        delete args_n;

        this->add(dPhi_n, 1.0, phi_np1, -1.0, phi_n, true);
        // Compute errors
        VectorXd errors = dPhi_n.getNorms();
#ifdef HAVE_MPI
        //distribute errors among all orbitals
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE, &errors(0), 1, MPI_DOUBLE, mpiCommOrb);
#endif

        phi_n.setErrors(errors);
        err_o = errors.maxCoeff();
        err_t = sqrt(errors.dot(errors));
        err_p = calcPropertyError();
        this->orbError.push_back(err_t);
        converged = checkConvergence(err_o, err_p);

        // Compute Fock matrix
        MatrixXd F_np1 = F_n + calcFockMatrixUpdate();
        phi_n.clear();
        dPhi_n.clear();
        fock.clear();

        // Rotate orbitals
        orthonormalize(fock, F_np1, phi_np1);

        // Update orbitals and Fock matrix
        int nOrbs = phi_np1.size();
        MatrixXd U = MatrixXd::Identity(nOrbs, nOrbs);
        F_n = F_np1;
        this->add.rotate(phi_n, U, phi_np1);
        phi_np1.clear();

        timer.stop();
        printOrbitals(F_n.diagonal(), phi_n, 0);
        printProperty();
        printTimer(timer.getWallTime());

        if (converged) break;
    }
    this->add.setPrecision(this->orbPrec[2]/10.0);
    if (this->canonical) {
        diagonalize(fock, F_n, phi_n);
    } else {
        localize(fock, F_n, phi_n);
    }

    printConvergence(converged);
    return converged;
}

MatrixXd EnergyOptimizer::calcFockMatrixUpdate() {
    if (this->fOper_np1 == 0) MSG_FATAL("Operator not initialized");

    OrbitalVector &phi_n = *this->orbitals_n;
    OrbitalVector &phi_np1 = *this->orbitals_np1;
    OrbitalVector &dPhi_n = *this->dOrbitals_n;

    TelePrompter::printHeader(0,"Computing Fock matrix update");

    Timer timer;
    IdentityOperator I;
    I.setup(getOrbitalPrecision());
    MatrixXd dS_1 = I(dPhi_n, phi_n);
    MatrixXd dS_2 = I(phi_np1, dPhi_n);
    I.clear();

    NuclearPotential *v_n = this->fOper_n->getNuclearPotential();
    CoulombOperator *j_n = this->fOper_n->getCoulombOperator();
    ExchangeOperator *k_n = this->fOper_n->getExchangeOperator();
    XCOperator *xc_n = this->fOper_n->getXCOperator();

    int Ni = phi_np1.size();
    int Nj = dPhi_n.size();
    MatrixXd dV_n = MatrixXd::Zero(Ni,Nj);
    {   // Nuclear potential matrix is computed explicitly
        Timer timer;

#ifdef HAVE_MPI

        OrbitalVector orbVecChunk_i(0); //to store adresses of own i_orbs
        OrbitalVector orbVecChunk_j(0); //to store adresses of own j_orbs
        OrbitalVector rcvOrbs(0);       //to store adresses of received orbitals

	vector<int> orbsIx;             //to store own orbital indices    
	int rcvOrbsIx[workOrbVecSize];  //to store received orbital indices

        //make vector with adresses of own orbitals
        for (int ix = mpiOrbRank; ix < Ni; ix += mpiOrbSize) {
            orbVecChunk_i.push_back(phi_np1.getOrbital(ix));//i orbitals
            orbsIx.push_back(ix);
        }
        for (int jx = mpiOrbRank; jx < Nj; jx += mpiOrbSize)
            orbVecChunk_j.push_back(dPhi_n.getOrbital(jx));//j orbitals

        for (int iter = 0; iter >= 0; iter++) {
            //get a new chunk from other processes
            orbVecChunk_i.getOrbVecChunk(orbsIx, rcvOrbs, rcvOrbsIx, Ni, iter);

            //Only one process does the computations. j orbitals always local
            MatrixXd resultChunk = MatrixXd::Zero(rcvOrbs.size(), orbVecChunk_j.size());
            resultChunk = (*v_n)(rcvOrbs, orbVecChunk_j);

            //copy results into final matrix
            int j = 0;
            for (int jx = mpiOrbRank;  jx < Nj; jx += mpiOrbSize) {
                for (int ix = 0; ix < rcvOrbs.size(); ix++) {
                    dV_n(rcvOrbsIx[ix],jx) += resultChunk(ix,j);
                }
                j++;
            }
            rcvOrbs.clearVec(false);//reset to zero size orbital vector
        }

        //clear orbital adresses, not the orbitals
        orbVecChunk_i.clearVec(false);
        orbVecChunk_j.clearVec(false);
        workOrbVec.clear();

        MPI_Allreduce(MPI_IN_PLACE, &dV_n(0,0), Ni*Nj,
                      MPI_DOUBLE, MPI_SUM, mpiCommOrb);

#else
        dV_n = (*v_n)(phi_np1, dPhi_n);
#endif

        timer.stop();
        double t = timer.getWallTime();
        TelePrompter::printDouble(0, "Nuclear potential matrix", t);
    }

    MatrixXd F_n;
    {   // Computing potential matrix excluding nuclear part
        Timer timer;
        FockOperator fock_n(0, 0, j_n, k_n, xc_n);
        F_n = fock_n(phi_np1, phi_n);
        timer.stop();
        double t = timer.getWallTime();
        TelePrompter::printDouble(0, "Fock matrix n", t);
    }

    {   // The n+1 Fock operator needs orthonormalized orbitals
        MatrixXd U = calcOrthonormalizationMatrix(phi_np1);
        this->add.rotate(phi_np1, U);
    }

    CoulombOperator *j_np1 = this->fOper_np1->getCoulombOperator();
    ExchangeOperator *k_np1 = this->fOper_np1->getExchangeOperator();
    XCOperator *xc_np1 = this->fOper_np1->getXCOperator();

    println(0,"                                                            ");
    // Do not setup exchange, it must be applied on the fly anyway
    if (j_np1 != 0) j_np1->setup(getOrbitalPrecision());
    if (k_np1 != 0) k_np1->ExchangeOperator::setup(getOrbitalPrecision());
    if (xc_np1 != 0) xc_np1->setup(getOrbitalPrecision());
    println(0,"                                                            ");

    MatrixXd F_np1;
    {   // Computing potential matrix excluding nuclear part
        Timer timer;
        FockOperator fock_np1(0, 0, j_np1, k_np1, xc_np1);
        MatrixXd F_1 = fock_np1(phi_n, phi_n);
        MatrixXd F_2 = fock_np1(phi_n, dPhi_n);
        fock_np1.clear();

        F_np1 = F_1 + F_2 + F_2.transpose();
        //MatrixXd F_3 = f_np1(*this->dPhi_n, *this->phi_n);
        //MatrixXd F_4 = f_np1(*this->dPhi_n, *this->dPhi_n);
        //MatrixXd F_np1 = F_1 + F_2 + F_3 + F_4;
        timer.stop();
        double t = timer.getWallTime();
        TelePrompter::printDouble(0, "Fock matrix n+1", t);
    }

    // Re-computing non-orthogonal phi_np1
    phi_np1.clear();
    this->add(phi_np1, 1.0, phi_n, 1.0, dPhi_n, true);

    MatrixXd L = this->helmholtz->getLambda().asDiagonal();
    MatrixXd dF_1 = dS_1*(*this->fMat_n);
    MatrixXd dF_2 = dS_2*L;
    MatrixXd dF_3 = F_np1 - F_n;

    // Adding up the pieces
    MatrixXd dF_n = dV_n + dF_1 + dF_2 + dF_3;

    // Symmetrizing
    MatrixXd sym = dF_n + dF_n.transpose();
    dF_n = 0.5 * sym;

    timer.stop();
    TelePrompter::printFooter(0, timer, 2);
    return dF_n;
}
