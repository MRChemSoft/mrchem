#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"

#include "OrbitalOptimizer.h"
#include "HelmholtzVector.h"
#include "KineticOperator.h"
#include "FockOperator.h"
#include "Accelerator.h"
#include "Orbital.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief constructor
 *
 * @param h: Helmholtz operators
 * @param k: Iterative accelerator
 *
 * SCF solver will NOT take ownership of the HelmholtzVector or the Accelerator,
 * so the original objects must be taken care of externally (do not delete until
 * SCF goes out of scope). Fock matrix, FockOperator and OrbitalVector are not
 * initialized at this stage, so the SCF solver needs to be "setup()" before
 * "optimize()".
 */
OrbitalOptimizer::OrbitalOptimizer(HelmholtzVector &h, Accelerator *k)
        : GroundStateSolver(h),
          kain(k) {
}

/** @brief Prepare solver for optimization
 *
 * @param fock: FockOperator defining the SCF equations
 * @param Phi: Orbitals to optimize
 * @param F: Fock matrix
 *
 * SCF solver will NOT take ownership of the input, so these objects must be taken
 * care of externally (do not delete until SCF goes out of scope).
 */
void OrbitalOptimizer::setup(FockOperator &fock,
                             OrbitalVector &Phi,
                             ComplexMatrix &F) {
    this->fMat_n = &F;
    this->fOper_n = &fock;
    this->orbitals_n = &Phi;
}

/** @brief Clear solver after optimization
 *
 * Clear pointers that was set during setup, and reset the precision parameter
 * (only the current precision orbPrec[0], not the boundary values orbPrec[1,2]).
 * Solver can be re-used after another setup.
 */
void OrbitalOptimizer::clear() {
    this->fMat_n = 0;
    this->fOper_n = 0;
    this->orbitals_n = 0;
    if (this->kain != 0) this->kain->clear();
    resetPrecision();
}

/** @brief Run orbital optimization
 *
 * Optimize orbitals until convergence thresholds are met. This algorithm computes
 * the Fock matrix explicitly using the kinetic energy operator, and uses a KAIN
 * accelerator to improve convergence. Diagonalization or localization may be performed
 * during the SCF iterations. Main points of the algorithm:
 *
 * Pre SCF: setup Fock operator and compute Fock matrix
 *
 *  1) Diagonalize/localize orbitals
 *  2) Compute current SCF energy
 *  3) Apply Helmholtz operator on all orbitals
 *  4) Orthonormalize orbitals (Löwdin)
 *  5) Compute orbital updates
 *  6) Compute KAIN update
 *  7) Compute errors and check for convergence
 *  8) Add orbital updates
 *  9) Orthonormalize orbitals (Löwdin)
 * 10) Setup Fock operator
 * 11) Compute Fock matrix
 *
 * Post SCF: diagonalize/localize orbitals
 */
bool OrbitalOptimizer::optimize() {
    ComplexMatrix &F = *this->fMat_n;
    FockOperator &fock = *this->fOper_n;
    OrbitalVector &Phi_n = *this->orbitals_n;
    HelmholtzVector &H = *this->helmholtz;

    double orb_prec = this->orbPrec[0];
    double err_o = orbital::get_errors(Phi_n).maxCoeff();
    double err_t = 1.0;
    double err_p = 1.0;

    fock.setup(orb_prec);
    F = fock.kinetic()(Phi_n, Phi_n) + fock.potential()(Phi_n, Phi_n);

    int nIter = 0;
    bool converged = false;
    while(nIter++ < this->maxIter or this->maxIter < 0) {
        // Initialize SCF cycle
        Timer timer;
        printCycleHeader(nIter);
        orb_prec = adjustPrecision(err_o);

        // Rotate orbitals
        if (needLocalization(nIter)) {
            ComplexMatrix U = orbital::localize(orb_prec, Phi_n);
            fock.rotate(U);
            F = U*F*U.adjoint();
            if (this->kain != 0) this->kain->clear();
        } else if (needDiagonalization(nIter)) {
            ComplexMatrix U = orbital::diagonalize(orb_prec, Phi_n, F);
            fock.rotate(U);
            if (this->kain != 0) this->kain->clear();
        }

        // Compute electronic energy
        double E = calcProperty();
        this->property.push_back(E);

        // Setup Helmholtz operators and argument
        H.setup(orb_prec, F.real().diagonal());
        ComplexMatrix L = H.getLambdaMatrix();
        OrbitalVector Psi_n = setupHelmholtzArguments(fock, L-F, Phi_n, true);

        // Apply Helmholtz operators
        OrbitalVector Phi_np1 = H(Psi_n);
        orbital::free(Psi_n);
        H.clear();

        ComplexMatrix U = orbital::orthonormalize(orb_prec, Phi_np1);
        F = U*F*U.adjoint();

        // Compute orbital updates
        OrbitalVector dPhi_n = orbital::add(1.0, Phi_np1, -1.0, Phi_n);
        orbital::free(Phi_np1);

        // Employ KAIN accelerator
        if (this->kain != 0) this->kain->accelerate(orb_prec, Phi_n, dPhi_n);

        // Compute errors
        DoubleVector errors = orbital::get_norms(dPhi_n);
        mpi::allreduce_vector(errors, mpi::comm_orb);

        err_o = errors.maxCoeff();
        err_t = errors.norm();
        err_p = calcPropertyError();
        this->orbError.push_back(err_t);
        converged = checkConvergence(err_o, err_p);

        // Update orbitals
        Phi_np1 = orbital::add(1.0, Phi_n, 1.0, dPhi_n);
        orbital::free(Phi_n);
        orbital::free(dPhi_n);
        Phi_n = Phi_np1;

        orbital::orthonormalize(orb_prec, Phi_n);
        orbital::set_errors(Phi_n, errors);

        // Compute Fock matrix
        fock.setup(orb_prec);
        F = fock.kinetic()(Phi_n, Phi_n) + fock.potential()(Phi_n, Phi_n);

        // Finalize SCF cycle
        timer.stop();
        printOrbitals(F.real().diagonal(), Phi_n, 0);
        printProperty();
        printCycleFooter(timer.getWallTime());

        if (converged) break;
    }

    if (this->kain != 0) this->kain->clear();
    fock.clear();

    if (this->canonical) {
        orbital::diagonalize(orb_prec/10, Phi_n, F);
    } else {
        ComplexMatrix U = orbital::localize(orb_prec/10, Phi_n);
        F = U*F*U.adjoint();
    }

    printConvergence(converged);
    return converged;
}

} //namespace mrchem
