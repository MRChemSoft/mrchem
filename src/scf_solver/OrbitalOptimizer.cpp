#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"

#include "OrbitalOptimizer.h"
#include "HelmholtzVector.h"
#include "FockOperator.h"
#include "Accelerator.h"
#include "Orbital.h"

using namespace std;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

OrbitalOptimizer::OrbitalOptimizer(HelmholtzVector &h, Accelerator *k)
    : GroundStateSolver(h),
      kain(k) {
}

OrbitalOptimizer::~OrbitalOptimizer() {
    this->kain = 0;
}

void OrbitalOptimizer::setup(FockOperator &fock,
                             OrbitalVector &Phi,
                             ComplexMatrix &F) {
    this->fMat_n = &F;
    this->fOper_n = &fock;

    this->orbitals_n = &Phi;
}

void OrbitalOptimizer::clear() {
    this->fMat_n = 0;
    this->fOper_n = 0;
    this->orbitals_n = 0;

    if (this->kain != 0) this->kain->clear();
    resetPrecision();
}

bool OrbitalOptimizer::optimize() {
    ComplexMatrix &F = *this->fMat_n;
    FockOperator &fock = *this->fOper_n;
    OrbitalVector &Phi_n = *this->orbitals_n;
    HelmholtzVector &H = *this->helmholtz;

    double orb_prec = getOrbitalPrecision();
    double err_o = orbital::get_errors(Phi_n).maxCoeff();
    double err_t = 1.0;
    double err_p = 1.0;

    fock.setup(orb_prec);
    F = fock(Phi_n, Phi_n);

    int nIter = 0;
    bool converged = false;
    while(nIter++ < this->maxIter or this->maxIter < 0) {
        // Initialize SCF cycle
        Timer timer;
        printCycle(nIter);
        adjustPrecision(err_o);
        orb_prec = getOrbitalPrecision();

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
        ComplexVector lambda = H.getLambda().cast<ComplexDouble>();
        ComplexMatrix L = lambda.asDiagonal();
        bool adjoint = false;
        bool clearFock = false;
        if (mpi::orb_size > 1) clearFock = true;
        OrbitalVector Psi_n = setupHelmholtzArguments(fock, L-F, Phi_n, adjoint, clearFock);

        // Apply Helmholtz operators
        OrbitalVector Phi_np1 = H(Psi_n);
        orbital::free(Psi_n);
        if (mpi::orb_size > 1) H.clear();

        if (not clearFock) fock.clear();
        ComplexMatrix U = orbital::orthonormalize(orb_prec, Phi_np1);
        F = U*F*U.adjoint();

        // Compute orbital updates
        OrbitalVector dPhi_n = orbital::add(1.0, Phi_np1, -1.0, Phi_n);
        orbital::free(Phi_np1);

        // Employ KAIN accelerator
        if (this->kain != 0) this->kain->accelerate(orb_prec, Phi_n, dPhi_n);

        // Compute errors
        int Ni = dPhi_n.size();
        DoubleVector errors = DoubleVector::Zero(Ni);
        for (int i = 0; i < Ni; i++) {
            if (mpi::my_orb(dPhi_n[i])) errors(i) = dPhi_n[i].norm();
        }

#ifdef HAVE_MPI
        //distribute errors among all orbitals
        MPI_Allreduce(MPI_IN_PLACE, &errors(0), Ni, MPI_DOUBLE, MPI_SUM, mpi::comm_orb);
#endif

        err_o = errors.maxCoeff();
        err_t = sqrt(errors.dot(errors));
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
        F = fock(Phi_n, Phi_n);

        // Finalize SCF cycle
        timer.stop();
        printOrbitals(F.real().diagonal(), Phi_n, 0);
        printProperty();
        printTimer(timer.getWallTime());

        if (converged) break;
    }
    if (this->kain != 0) this->kain->clear();
    fock.clear();

    if (this->canonical) {
        orbital::diagonalize(orb_prec, Phi_n, F);
    } else {
        ComplexMatrix U = orbital::localize(orb_prec, Phi_n);
        F = U*F*U.adjoint();
    }

    printConvergence(converged);
    return converged;
}

/** Prints the number of trees and nodes kept in the solver at the given moment */
void OrbitalOptimizer::printTreeSizes() const {
    NOT_IMPLEMENTED_ABORT;
    /*
    void printTreeSizes() const;
    TelePrompter::printHeader(0, "Printing Tree Sizes");

    int nNodes = 0;
    if (this->fOper_n != 0) nNodes += this->fOper_n->printTreeSizes();
    if (this->orbitals_n != 0) nNodes += this->orbitals_n->printTreeSizes();
    if (this->kain != 0) nNodes += this->kain->printTreeSizes();

    TelePrompter::printSeparator(0, '-');
    println(0," Total number of nodes                   " << setw(18) << nNodes);
    TelePrompter::printSeparator(0, '=', 2);
*/
}

} //namespace mrchem
