#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"

#include "EnergyOptimizer.h"
#include "HelmholtzVector.h"
#include "Orbital.h"
#include "orbital_utils.h"

#include "FockOperator.h"
#include "NuclearOperator.h"
#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "XCOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief constructor
 *
 * @param h: Helmholtz operators
 *
 * SCF solver will NOT take ownership of the HelmholtzVector, so the original object
 * must be taken care of externally (do not delete until SCF goes out of scope).
 * Fock matrix, FockOperator and OrbitalVector are not initialized at this stage,
 * so the SCF solver needs to be "setup()" before "optimize()".
 */
EnergyOptimizer::EnergyOptimizer(HelmholtzVector &h)
        : GroundStateSolver(h),
          fOper_np1(0),
          orbitals_np1(0) {
}

/** @brief Prepare solver for optimization
 *
 * @param fock: FockOperator defining the SCF equations
 * @param Phi: Orbitals to optimize
 * @param F: Fock matrix
 * @param fock_np1: Next iteration FockOperator
 * @param Phi_np1: Next iteration orbitals
 *
 * SCF solver will NOT take ownership of the input, so these objects must be taken
 * care of externally (do not delete until SCF goes out of scope).
 */
void EnergyOptimizer::setup(FockOperator &fock, OrbitalVector &Phi, ComplexMatrix &F,
                            FockOperator &fock_np1, OrbitalVector &Phi_np1) {
    this->fMat_n = &F;
    this->fOper_n = &fock;
    this->fOper_np1 = &fock_np1;
    this->orbitals_n = &Phi;
    this->orbitals_np1 = &Phi_np1;
}

/** @brief Clear solver after optimization
 *
 * Clear pointers that was set during setup, and reset the precision parameter
 * (only the current precision orbPrec[0], not the bounaary values orbPrec[1,2]).
 * Solver can be re-used after another setup.
 */
void EnergyOptimizer::clear() {
    this->fMat_n = 0;
    this->fOper_n = 0;
    this->fOper_np1 = 0;
    this->orbitals_n = 0;
    this->orbitals_np1 = 0;
    resetPrecision();
}

/** @brief Run energy optimization
 *
 * Optimize orbitals until convergence thresholds are met. This algorithm does NOT
 * use any kinetic energy operator, but computes the Fock matrix as a potential
 * update from the previous iteration. No KAIN accelerator and no diagonalization
 * or localization within the SCF iterations. Main points of the algorithm:
 *
 * Pre SCF: diagonalize/localize orbitals
 *
 * 1) Setup Fock operator
 * 2) Compute current SCF energy
 * 3) Apply Helmholtz operator on all orbitals
 * 4) Compute orbital updates
 * 5) Compute errors and check for convergence
 * 6) Compute Fock matrix update
 * 7) Orthonormalize orbitals (Löwdin)
 *
 * Post SCF: diagonalize/localize orbitals
 */
bool EnergyOptimizer::optimize() {
    ComplexMatrix &F_n = *this->fMat_n;
    FockOperator &fock = *this->fOper_n;
    OrbitalVector &Phi_n = *this->orbitals_n;
    OrbitalVector &Phi_np1 = *this->orbitals_np1;
    HelmholtzVector &H = *this->helmholtz;

    double orb_prec = this->orbPrec[0];
    double err_o = orbital::get_errors(Phi_n).maxCoeff();
    double err_t = 1.0;
    double err_p = 1.0;

    if (this->canonical) {
        orbital::diagonalize(orb_prec, Phi_n, F_n);
    } else {
        ComplexMatrix U = orbital::localize(orb_prec, Phi_n);
        F_n = U*F_n*U.adjoint();
    }

    int nIter = 0;
    bool converged = false;
    while(nIter++ < this->maxIter or this->maxIter < 0) {
        // Initialize SCF cycle
        Timer timer;
        printCycleHeader(nIter);
        orb_prec = adjustPrecision(err_o);

        // Compute electronic energy
        fock.setup(orb_prec);
        double E = calcProperty();
        this->property.push_back(E);

        // Setup Helmholtz operators and argument
        H.setup(orb_prec, F_n.real().diagonal());
        ComplexMatrix L_n = H.getLambdaMatrix();
        OrbitalVector Psi_n = setupHelmholtzArguments(fock, L_n-F_n, Phi_n, false);

        // Apply Helmholtz operators
        Phi_np1 = H(Psi_n);
        orbital::free(Psi_n);

        // Compute orbital updates
        OrbitalVector dPhi_n = orbital::add(1.0, Phi_np1, -1.0, Phi_n);

        // Compute orbital errors
        DoubleVector errors = orbital::get_norms(dPhi_n);
        mpi::allreduce_vector(errors, mpi::comm_orb);

        err_o = errors.maxCoeff();
        err_t = errors.norm();
        err_p = calcPropertyError();
        this->orbError.push_back(err_t);
        converged = checkConvergence(err_o, err_p);

        // Compute Fock matrix
        ComplexMatrix F_np1 = F_n + calcFockMatrixUpdate(orb_prec, dPhi_n);
        H.clear();
        orbital::free(Phi_n);
        orbital::free(dPhi_n);
        fock.clear();

        // Rotate orbitals
        ComplexMatrix U = orbital::calc_lowdin_matrix(Phi_np1);
        Phi_n = orbital::linear_combination(U, Phi_np1, orb_prec);
        F_n = U*F_np1*U.adjoint();
        orbital::free(Phi_np1);
        orbital::set_errors(Phi_n, errors);

        timer.stop();
        printOrbitals(F_n.real().diagonal(), Phi_n, 0);
        printProperty();
        printCycleFooter(timer.getWallTime());

        if (converged) break;
    }

    if (this->canonical) {
        orbital::diagonalize(orb_prec/10, Phi_n, F_n);
    } else {
        ComplexMatrix U = orbital::localize(orb_prec/10, Phi_n);
        F_n = U*F_n*U.adjoint();
    }

    printConvergence(converged);
    return converged;
}

/** @brief Compute Fock matrix update
 *
 * The Fock matrix is computed as a direct update from the previous iteration.
 * This update is exact, provided that the orbital relation between the orbitals
 * at iterations n and n+1 is exactly the application of the Helmholtz operator.
 *
 * \Delta F = \Delta S_1 F + \Delta S_2 \Lambda + \Delta F_{pot}
 *
 * (\Delta S_1^n)_{ij} = <\Delta \phi_i^{n}   | \phi_j^n>
 * (\Delta S_2^n)_{ij} = <\Delta \phi_i^{n+1} | \Delta \phi_j^n>
 *
 * (\Delta F_^n{pot})_{ij} = <\phi_i^{n+1} | V^n | \Delta \phi_j^n>
 *                         + <\phi_i^{n+1} | \Delta V^n | \phi_j^n>
 *
 */
ComplexMatrix EnergyOptimizer::calcFockMatrixUpdate(double prec, OrbitalVector &dPhi_n) {
    if (this->fOper_np1 == 0) MSG_FATAL("Operator not initialized");

    OrbitalVector &Phi_n = *this->orbitals_n;
    OrbitalVector &Phi_np1 = *this->orbitals_np1;

    Printer::printHeader(0,"Computing Fock matrix update");

    Timer timer;
    ComplexMatrix dS_1 = orbital::calc_overlap_matrix(dPhi_n, Phi_n);
    ComplexMatrix dS_2 = orbital::calc_overlap_matrix(Phi_np1, dPhi_n);

    NuclearOperator  *v_n = this->fOper_n->getNuclearOperator();
    CoulombOperator  *j_n = this->fOper_n->getCoulombOperator();
    ExchangeOperator *k_n = this->fOper_n->getExchangeOperator();
    XCOperator      *xc_n = this->fOper_n->getXCOperator();

        /*
    int Ni = Phi_np1.size();
    int Nj = dPhi_n.size();
    ComplexMatrix dV_n = ComplexMatrix::Zero(Ni,Nj);
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
            ComplexMatrix resultChunk = ComplexMatrix::Zero(rcvOrbs.size(), orbVecChunk_j.size());
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
        Printer::printDouble(0, "Nuclear potential matrix", t, 5);
    }
    */

    ComplexMatrix dV_n;
    {   // Nuclear potential matrix is computed explicitly
        Timer timer;
        dV_n = (*v_n)(Phi_np1, dPhi_n);
        timer.stop();
        double t = timer.getWallTime();
        Printer::printDouble(0, "Nuclear potential matrix", t, 5);
    }

    ComplexMatrix F_n;
    {   // Computing potential matrix excluding nuclear part
        Timer timer;
        FockOperator fock_n(0, 0, j_n, k_n, xc_n);
        F_n = fock_n(Phi_np1, Phi_n);
        timer.stop();
        double t = timer.getWallTime();
        Printer::printDouble(0, "Fock matrix n", t, 5);
    }

    {   // The n+1 Fock operator needs orthonormalized orbitals
        orbital::orthonormalize(prec, Phi_np1);
    }

    CoulombOperator  *j_np1 = this->fOper_np1->getCoulombOperator();
    ExchangeOperator *k_np1 = this->fOper_np1->getExchangeOperator();
    XCOperator      *xc_np1 = this->fOper_np1->getXCOperator();

    println(0,"                                                            ");
    // Do not setup internal exchange, it must be applied on the fly anyway
    if (j_np1 != 0)   j_np1->setup(prec);
    if (k_np1 != 0)   k_np1->setup(prec);
    if (xc_np1 != 0) xc_np1->setup(prec);
    println(0,"                                                            ");

    ComplexMatrix F_np1;
    {   // Computing potential matrix excluding nuclear part
        Timer timer;
        FockOperator fock_np1(0, 0, j_np1, k_np1, xc_np1);
        ComplexMatrix F_1 = fock_np1(Phi_n, Phi_n);
        ComplexMatrix F_2 = fock_np1(Phi_n, dPhi_n);
        fock_np1.clear();

        F_np1 = F_1 + F_2 + F_2.transpose();
        //ComplexMatrix F_3 = f_np1(*this->dPhi_n, *this->phi_n);
        //ComplexMatrix F_4 = f_np1(*this->dPhi_n, *this->dPhi_n);
        //ComplexMatrix F_np1 = F_1 + F_2 + F_3 + F_4;
        timer.stop();
        double t = timer.getWallTime();
        Printer::printDouble(0, "Fock matrix n+1", t, 5);
    }

    // Re-computing non-orthogonal phi_np1
    orbital::free(Phi_np1);
    Phi_np1 = orbital::add(1.0, Phi_n, 1.0, dPhi_n);

    ComplexMatrix L = this->helmholtz->getLambdaMatrix();
    ComplexMatrix dF_1 = dS_1*(*this->fMat_n);
    ComplexMatrix dF_2 = dS_2*L;
    ComplexMatrix dF_3 = F_np1 - F_n;

    // Adding up the pieces
    ComplexMatrix dF_n = dV_n + dF_1 + dF_2 + dF_3;

    // Symmetrizing
    ComplexMatrix sym = dF_n + dF_n.transpose();
    dF_n = 0.5 * sym;

    timer.stop();
    Printer::printFooter(0, timer, 2);
    return dF_n;
}

} //namespace mrchem
