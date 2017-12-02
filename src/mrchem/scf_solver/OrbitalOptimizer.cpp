#include "OrbitalOptimizer.h"
#include "OrbitalVector.h"
#include "FockOperator.h"
#include "HelmholtzOperatorSet.h"
#include "Accelerator.h"

using namespace std;
using namespace Eigen;

OrbitalOptimizer::OrbitalOptimizer(HelmholtzOperatorSet &h,
                                   Accelerator *k)
    : GroundStateSolver(h),
      kain(k) {
}

OrbitalOptimizer::~OrbitalOptimizer() {
    this->kain = 0;
}

void OrbitalOptimizer::setup(FockOperator &fock,
                             OrbitalVector &phi,
                             MatrixXd &F) {
    this->fMat_n = &F;
    this->fOper_n = &fock;

    this->orbitals_n = &phi;
    this->orbitals_np1 = new OrbitalVector(phi);
    this->dOrbitals_n = new OrbitalVector(phi);
}

void OrbitalOptimizer::clear() {
    this->fMat_n = 0;
    this->fOper_n = 0;

    delete this->orbitals_np1;
    delete this->dOrbitals_n;

    this->orbitals_n = 0;
    this->orbitals_np1 = 0;
    this->dOrbitals_n = 0;

    if (this->kain != 0) this->kain->clear();
    resetPrecision();
}

bool OrbitalOptimizer::optimize() {
    MatrixXd &F = *this->fMat_n;
    FockOperator &fock = *this->fOper_n;
    OrbitalVector &phi_n = *this->orbitals_n;
    OrbitalVector &phi_np1 = *this->orbitals_np1;
    OrbitalVector &dPhi_n = *this->dOrbitals_n;
    HelmholtzOperatorSet &H = *this->helmholtz;

    double orb_prec = getOrbitalPrecision();
    double err_o = phi_n.getErrors().maxCoeff();
    double err_t = 1.0;
    double err_p = 1.0;

    fock.setup(orb_prec);
    F = fock(phi_n, phi_n);

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
            localize(fock, F, phi_n);
            if (this->kain != 0) this->kain->clear();
        } else if (needDiagonalization(nIter)) {
            diagonalize(fock, F, phi_n);
            if (this->kain != 0) this->kain->clear();
        }

        // Compute electronic energy
        double E = calcProperty();
        this->property.push_back(E);

        // Setup Helmholtz operators and argument
        H.setup(orb_prec, F.diagonal());
        MatrixXd L = H.getLambda().asDiagonal();
	bool adjoint = false;
	bool clearFock = false;
	if (mpiOrbSize > 1) clearFock = true;
        OrbitalVector *args_n = setupHelmholtzArguments(fock, L-F, phi_n, adjoint, clearFock);

        // Apply Helmholtz operators
        H(phi_np1, *args_n);
        delete args_n;
        if (mpiOrbSize > 1) H.clear();

        if (not clearFock) fock.clear();
        orthonormalize(fock, F, phi_np1);

        // Compute orbital updates
        this->add(dPhi_n, 1.0, phi_np1, -1.0, phi_n, true);
        phi_np1.clear();

        // Employ KAIN accelerator
        if (this->kain != 0) this->kain->accelerate(orb_prec, phi_n, dPhi_n);

        // Compute errors
        int Ni = dPhi_n.size();
        VectorXd errors = VectorXd::Zero(Ni);
        for (int i = mpiOrbRank; i < Ni; i += mpiOrbSize){
            errors(i) = sqrt(dPhi_n.getOrbital(i).getSquareNorm());
        }

#ifdef HAVE_MPI
        //distribute errors among all orbitals
        MPI_Allreduce(MPI_IN_PLACE, &errors(0), Ni, MPI_DOUBLE, MPI_SUM, mpiCommOrb);
#endif

        err_o = errors.maxCoeff();
        err_t = sqrt(errors.dot(errors));
        err_p = calcPropertyError();
        this->orbError.push_back(err_t);
        converged = checkConvergence(err_o, err_p);

        // Update orbitals
        this->add.inPlace(phi_n, 1.0, dPhi_n);
        dPhi_n.clear();

        orthonormalize(fock, F, phi_n);
        phi_n.setErrors(errors);

        // Compute Fock matrix
        fock.setup(orb_prec);
        F = fock(phi_n, phi_n);

        // Finalize SCF cycle
        timer.stop();
        printOrbitals(F.diagonal(), phi_n, 0);
        printProperty();
        printTimer(timer.getWallTime());

        if (converged) break;
    }
    if (this->kain != 0) this->kain->clear();
    fock.clear();

    this->add.setPrecision(this->orbPrec[2]/10.0);
    if (this->canonical) {
        diagonalize(fock, F, phi_n);
    } else {
        localize(fock, F, phi_n);
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
    if (this->orbitals_np1 != 0) nNodes += this->orbitals_np1->printTreeSizes();
    if (this->dOrbitals_n != 0) nNodes += this->dOrbitals_n->printTreeSizes();
    if (this->kain != 0) nNodes += this->kain->printTreeSizes();

    TelePrompter::printSeparator(0, '-');
    println(0," Total number of nodes                   " << setw(18) << nNodes);
    TelePrompter::printSeparator(0, '=', 2);
*/
}
