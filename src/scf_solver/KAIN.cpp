#include "KAIN.h"
#include "Orbital.h"
#include "OrbitalVector.h"
#include "eigen_disable_warnings.h"

using namespace std;
using namespace Eigen;

/** Calculates the A matrices and b vectors for each orbital individually.
 *
  * \f$ A_{ij} = \langle x^n - x^i | f(x^n) - f(x^j) \rangle \f$
  * \f$ b_{i} = \langle x^n - x^i | f(x^n) \rangle \f$
  *
  * If Fock matrix is included this is treated as an additional ''orbital''
  * and the return vectors have size nOrbs + 1. Frobenius inner product
  * used for the Fock matrix. If separateOrbitals is false the A's and b's
  * are later collected to single entities. */
void KAIN::setupLinearSystem() {
    Timer timer;
    int nHistory = this->orbitals.size() - 1;
    if (nHistory < 1) {
        MSG_FATAL("Not enough history to setup system of equations");
    }

    vector<MatrixXd *> A_matrices;
    vector<VectorXd *> b_vectors;

    int nOrbitals = this->orbitals[nHistory]->size();
    for (int n = 0; n < nOrbitals; n++) {
        MatrixXd *orbA = new MatrixXd;
        VectorXd *orbB = new VectorXd;
        *orbA = MatrixXd::Zero(nHistory, nHistory);
        *orbB = VectorXd::Zero(nHistory);

	if (mpiOrbRank == n%mpiOrbSize) {
	    Orbital &phi_m = this->orbitals[nHistory]->getOrbital(n);
	    Orbital &fPhi_m = this->dOrbitals[nHistory]->getOrbital(n);

	    for (int i = 0; i < nHistory; i++) {
		Orbital &phi_i = this->orbitals[i]->getOrbital(n);
		Orbital dPhi_im(phi_i);
		this->add(dPhi_im, 1.0, phi_i, -1.0, phi_m, true);

		for (int j = 0; j < nHistory; j++) {
		    Orbital &fPhi_j = this->dOrbitals[j]->getOrbital(n);
		    Orbital dfPhi_jm(fPhi_j);
		    this->add(dfPhi_jm, 1.0, fPhi_j, -1.0, fPhi_m, true);

		    // Ref. Harrisons KAIN paper the following has the wrong sign,
		    // but we define the updates (lowercase f) with opposite sign.
		    complex<double> inner_prod = dPhi_im.dot(dfPhi_jm);
		    if (inner_prod.imag() > MachineZero) NOT_IMPLEMENTED_ABORT;
		    (*orbA)(i,j) -= inner_prod.real();
		}
		complex<double> inner_prod = dPhi_im.dot(fPhi_m);
		if (inner_prod.imag() > MachineZero) NOT_IMPLEMENTED_ABORT;
		(*orbB)(i) += inner_prod.real();
	    }
	}
        A_matrices.push_back(orbA);
        b_vectors.push_back(orbB);
    }

#ifdef HAVE_MPI
    //combine results from all processes
    //make a super matrix, so that all coefficients are consecutive
    MatrixXd orbsAB(nHistory, (nHistory+1)*nOrbitals);
    
    int jn = 0;
    for (int n = 0; n < nOrbitals; n++) {
	MatrixXd *orbA = A_matrices[n];
	for (int j = 0; j < nHistory; j++) {
	    for (int i = 0; i < nHistory; i++) {
		orbsAB(i,jn) = (*orbA)(i,j);
	    }
	    jn++;
	}
    }
    for (int n = 0; n < nOrbitals; n++) {
	VectorXd *orbB = b_vectors[n];
	for (int i = 0; i < nHistory; i++) {
	    orbsAB(i,jn) = (*orbB)(i);
	}
	jn++;
    }

    MPI_Allreduce(MPI_IN_PLACE, &orbsAB(0,0), (nHistory+1)*nHistory*nOrbitals,
		  MPI_DOUBLE, MPI_SUM, mpiCommOrb);
    
    jn = 0;
    for (int n = 0; n < nOrbitals; n++) {
	MatrixXd *orbA = A_matrices[n];
	for (int j = 0; j < nHistory; j++) {
	    for (int i = 0; i < nHistory; i++) {
		(*orbA)(i,j) = orbsAB(i,jn);
	    }
	    jn++;
	}
    }
    for (int n = 0; n < nOrbitals; n++) {
	VectorXd *orbB = b_vectors[n];
	for (int i = 0; i < nHistory; i++) {
	    (*orbB)(i) = orbsAB(i,jn);
	}
	jn++;
    }

#endif

//    // Fock matrix is treated as a whole using the Frobenius inner product
    if (this->orbitals.size() == this->fock.size()) {
        const MatrixXd &X_m = this->fock[nHistory];
        const MatrixXd &fX_m = this->dFock[nHistory];

        MatrixXd *fockA = new MatrixXd;
        *fockA = MatrixXd::Zero(nHistory, nHistory);

        VectorXd *fockB = new VectorXd;
        *fockB = VectorXd::Zero(nHistory);

        for (int i = 0; i < nHistory; i++) {
            const MatrixXd &X_i = this->fock[i];
            MatrixXd dX_im = X_i - X_m;
            for (int j = 0; j < nHistory; j++) {
                const MatrixXd &fX_j = this->dFock[j];
                MatrixXd dfX_jm = fX_j - fX_m;
                MatrixXd prod = dX_im.transpose()*dfX_jm;
                (*fockA)(i,j) -= prod.trace();
            }
            MatrixXd prod = dX_im.transpose()*fX_m;
            (*fockB)(i) += prod.trace();
        }
        A_matrices.push_back(fockA);
        b_vectors.push_back(fockB);
    }

    sortLinearSystem(A_matrices, b_vectors);

    for (int n = 0; n < A_matrices.size(); n++) {
        delete A_matrices[n];
        delete b_vectors[n];
    }
    timer.stop();
    double t = timer.getWallTime();
    TelePrompter::printDouble(0, "Setup linear system", t);
}

/** Compute the next step for orbitals and orbital updates given the
  * solution \f$ c \f$ of the linear problem \f$ Ac = b \f$.
  *
  * The next step \f$ \delta x^n \f$ is constructed as
  * \f$ \delta x^n = f(x^n) + \sum_{j=1}^m c_j[(x^j-x^n)+(f(x^j)-f(x^n))]\f$
  */
void KAIN::expandSolution(OrbitalVector &phi,
                          OrbitalVector &dPhi,
                          MatrixXd *F,
                          MatrixXd *dF) {
    Timer timer;
    int nHistory = this->orbitals.size() - 1;
    int nOrbitals = this->orbitals[nHistory]->size();

    phi.clear();
    dPhi.clear();

    // Orbitals are unchanged
    MatrixXd I = MatrixXd::Identity(phi.size(), phi.size());
    this->add.rotate(phi, I, *this->orbitals[nHistory]);

    int m = 0;
    for (int n = 0; n < nOrbitals; n++) {
        if (this->sepOrbitals) {
            m = n;
        }
	if (mpiOrbRank == n%mpiOrbSize) {
	    vector<complex<double> > totCoefs;
	    vector<Orbital *> totOrbs;

	    Orbital &phi_m = this->orbitals[nHistory]->getOrbital(n);
	    Orbital &fPhi_m = this->dOrbitals[nHistory]->getOrbital(n);
	    totCoefs.push_back(1.0);
	    totOrbs.push_back(&fPhi_m);

	    // Ref. Harrisons KAIN paper the following has the wrong sign,
	    // but we define the updates (lowercase f) with opposite sign
	    // (but not the orbitals themselves).
	    for (int j = 0; j < nHistory; j++) {
		vector<complex<double> > partCoefs;
		vector<Orbital *> partOrbs;

		Orbital &phi_j = this->orbitals[j]->getOrbital(n);
		partCoefs.push_back(1.0);
		partOrbs.push_back(&phi_j);

		Orbital &fPhi_j = this->dOrbitals[j]->getOrbital(n);
		partCoefs.push_back(1.0);
		partOrbs.push_back(&fPhi_j);

		partCoefs.push_back(-1.0);
		partOrbs.push_back(&phi_m);
		partCoefs.push_back(-1.0);
		partOrbs.push_back(&fPhi_m);

		Orbital *partStep = new Orbital(fPhi_m);
		this->add(*partStep, partCoefs, partOrbs, true);

		double c_j = (*this->c[m])(j);
		totCoefs.push_back(c_j);
		totOrbs.push_back(partStep);
	    }

	    Orbital &dPhi_n = dPhi.getOrbital(n);
	    this->add(dPhi_n, totCoefs, totOrbs, true);
	    for (int k = 0; k < totOrbs.size(); k++) {
		// First entry is the last orbital update and should not be deallocated,
		// all other entries are locally allocated partSteps that must be deleted
		if (k != 0) delete totOrbs[k];
		totOrbs[k] = 0;
	    }
	}
    }

    // Treat Fock matrix as a whole using Frobenius inner product
    if (this->fock.size() == this->orbitals.size()) {
        if (F == 0 or dF == 0) MSG_ERROR("Invalid fock matrix");

        MatrixXd X_m = this->fock[nHistory];
        const MatrixXd &fX_m = this->dFock[nHistory];
        MatrixXd fockStep = MatrixXd::Zero(nOrbitals, nOrbitals);
        fockStep = fX_m;
        m = this->c.size();
        for (int j = 0; j < nHistory; j++) {
            const MatrixXd &X_j = this->fock[j];
            const MatrixXd &fX_j = this->dFock[j];
            MatrixXd tmpX = MatrixXd::Zero(nOrbitals,nOrbitals);
            tmpX += X_j;
            tmpX -= X_m;
            tmpX += fX_j;
            tmpX -= fX_m;
            tmpX *= (*this->c[m-1])(j);
            fockStep += tmpX;
        }
        *F = X_m;
        *dF = fockStep;
    }
    timer.stop();
    double t = timer.getWallTime();
    TelePrompter::printDouble(0, "Expand solution", t);
}
