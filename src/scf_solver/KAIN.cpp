#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"

#include "KAIN.h"
#include "Orbital.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Calculates the A matrices and b vectors for each orbital individually.
 *
 * \f$ A_{ij} = \langle x^n - x^i | f(x^n) - f(x^j) \rangle \f$
 * \f$ b_{i} = \langle x^n - x^i | f(x^n) \rangle \f$
 *
 * If Fock matrix is included this is treated as an additional ''orbital''
 * and the return vectors have size nOrbs + 1. Frobenius inner product
 * used for the Fock matrix. If separateOrbitals is false the A's and b's
 * are later collected to single entities.
 */
void KAIN::setupLinearSystem() {
    Timer timer;
    int nHistory = this->orbitals.size() - 1;
    if (nHistory < 1) MSG_FATAL("Not enough history to setup system of equations");

    std::vector<DoubleMatrix> A_matrices;
    std::vector<DoubleVector> b_vectors;

    int nOrbitals = this->orbitals[nHistory].size();
    for (int n = 0; n < nOrbitals; n++) {
        DoubleMatrix orbA = DoubleMatrix::Zero(nHistory, nHistory);
        DoubleVector orbB = DoubleVector::Zero(nHistory);

        if (mpi::orb_rank == n%mpi::orb_size) {
            Orbital &phi_m = this->orbitals[nHistory][n];
            Orbital &fPhi_m = this->dOrbitals[nHistory][n];

            for (int i = 0; i < nHistory; i++) {
                Orbital &phi_i = this->orbitals[i][n];
                Orbital dPhi_im = orbital::add(1.0, phi_i, -1.0, phi_m);

                for (int j = 0; j < nHistory; j++) {
                    Orbital &fPhi_j = this->dOrbitals[j][n];
                    Orbital dfPhi_jm = orbital::add(1.0, fPhi_j, -1.0, fPhi_m);

                    // Ref. Harrisons KAIN paper the following has the wrong sign,
                    // but we define the updates (lowercase f) with opposite sign.
                    orbA(i,j) -= orbital::dot(dPhi_im, dfPhi_jm).real();
                }
                orbB(i) += orbital::dot(dPhi_im, fPhi_m).real();
            }
        }
        A_matrices.push_back(orbA);
        b_vectors.push_back(orbB);
    }

    /*
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
*/

    // Fock matrix is treated as a whole using the Frobenius inner product
    if (this->orbitals.size() == this->fock.size()) {
        const ComplexMatrix &X_m = this->fock[nHistory];
        const ComplexMatrix &fX_m = this->dFock[nHistory];

        ComplexMatrix fockA = ComplexMatrix::Zero(nHistory, nHistory);

        ComplexVector fockB = ComplexVector::Zero(nHistory);

        for (int i = 0; i < nHistory; i++) {
            const ComplexMatrix &X_i = this->fock[i];
            ComplexMatrix dX_im = X_i - X_m;
            for (int j = 0; j < nHistory; j++) {
                const ComplexMatrix &fX_j = this->dFock[j];
                ComplexMatrix dfX_jm = fX_j - fX_m;
                ComplexMatrix prod = dX_im.transpose()*dfX_jm;
                fockA(i,j) -= prod.trace();
            }
            ComplexMatrix prod = dX_im.transpose()*fX_m;
            fockB(i) += prod.trace();
        }
        A_matrices.push_back(fockA.real());
        b_vectors.push_back(fockB.real());
    }

    sortLinearSystem(A_matrices, b_vectors);

    timer.stop();
    double t = timer.getWallTime();
    Printer::printDouble(0, "Setup linear system", t, 5);
}

/** @brief Compute the next step for orbitals and orbital updates
 *
 * The next step \f$ \delta x^n \f$ is constructed from the solution
 * solution \f$ c \f$ of the linear problem \f$ Ac = b \f$ as:
 *
 * \f$ \delta x^n = f(x^n) + \sum_{j=1}^m c_j[(x^j-x^n)+(f(x^j)-f(x^n))]\f$
 */
void KAIN::expandSolution(double prec,
                          OrbitalVector &Phi,
                          OrbitalVector &dPhi,
                          ComplexMatrix *F,
                          ComplexMatrix *dF) {
    Timer timer;
    int nHistory = this->orbitals.size() - 1;
    int nOrbitals = this->orbitals[nHistory].size();

    Phi.clear();
    dPhi.clear();

    // Orbitals are unchanged
    Phi = orbital::deep_copy(this->orbitals[nHistory]);

    int m = 0;
    for (int n = 0; n < nOrbitals; n++) {
        if (this->sepOrbitals) {
            m = n;
        }
        if (mpi::orb_rank == n%mpi::orb_size) {
            std::vector<ComplexDouble> totCoefs;
            OrbitalVector totOrbs;

            Orbital &phi_m = this->orbitals[nHistory][n];
            Orbital &fPhi_m = this->dOrbitals[nHistory][n];
            totCoefs.push_back(1.0);
            totOrbs.push_back(fPhi_m);

            // Ref. Harrisons KAIN paper the following has the wrong sign,
            // but we define the updates (lowercase f) with opposite sign
            // (but not the orbitals themselves).
            for (int j = 0; j < nHistory; j++) {
                ComplexVector partCoefs(4);
                OrbitalVector partOrbs;

                partCoefs(0) = 1.0;
                Orbital &phi_j = this->orbitals[j][n];
                partOrbs.push_back(phi_j);

                partCoefs(1) = 1.0;
                Orbital &fPhi_j = this->dOrbitals[j][n];
                partOrbs.push_back(fPhi_j);

                partCoefs(2) = -1.0;
                partOrbs.push_back(phi_m);

                partCoefs(3) = -1.0;
                partOrbs.push_back(fPhi_m);

                Orbital partStep = orbital::multiply(partCoefs, partOrbs, prec);

                ComplexDouble c_j = this->c[m](j);
                totCoefs.push_back(c_j);
                totOrbs.push_back(partStep);
            }

            // std::vector -> ComplexVector
            ComplexVector coefsVec(totCoefs.size());
            for (int i = 0; i < totCoefs.size(); i++) coefsVec(i) = totCoefs[i];

            Orbital dPhi_n = orbital::multiply(coefsVec, totOrbs, prec);
            dPhi.push_back(dPhi_n);

            // First entry is the last orbital update and should not be deallocated,
            // all other entries are locally allocated partSteps that must be deleted
            totOrbs[0].clear();
            orbital::free(totOrbs);
        }
    }

    // Treat Fock matrix as a whole using Frobenius inner product
    if (this->fock.size() == this->orbitals.size()) {
        if (F == 0 or dF == 0) MSG_ERROR("Invalid fock matrix");

        ComplexMatrix X_m = this->fock[nHistory];
        const ComplexMatrix &fX_m = this->dFock[nHistory];
        ComplexMatrix fockStep = ComplexMatrix::Zero(nOrbitals, nOrbitals);
        fockStep = fX_m;
        m = this->c.size();
        for (int j = 0; j < nHistory; j++) {
            const ComplexMatrix &X_j = this->fock[j];
            const ComplexMatrix &fX_j = this->dFock[j];
            ComplexMatrix tmpX = ComplexMatrix::Zero(nOrbitals,nOrbitals);
            tmpX += X_j;
            tmpX -= X_m;
            tmpX += fX_j;
            tmpX -= fX_m;
            tmpX *= this->c[m-1](j);
            fockStep += tmpX;
        }
        *F = X_m;
        *dF = fockStep;
    }
    timer.stop();
    double t = timer.getWallTime();
    Printer::printDouble(0, "Expand solution", t, 5);
}

} //namespace mrchem
