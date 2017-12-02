#include "SCF.h"
#include "HelmholtzOperatorSet.h"
#include "FockOperator.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "Timer.h"

using namespace std;
using namespace Eigen;

extern MultiResolutionAnalysis<3> *MRA; // Global MRA
extern OrbitalVector workOrbVec;

SCF::SCF(HelmholtzOperatorSet &h)
    : maxIter(-1),
      rotation(0),
      canonical(true),
      orbThrs(-1.0),
      propThrs(-1.0),
      add(-1.0, MRA->getMaxScale()),
      helmholtz(&h) {
    this->orbPrec[0] = -1.0;
    this->orbPrec[1] = -1.0;
    this->orbPrec[2] = -1.0;
}

SCF::~SCF() {
    this->helmholtz = 0;
}

void SCF::setThreshold(double orb_thrs, double prop_thrs) {
    this->orbThrs = orb_thrs;
    this->propThrs = prop_thrs;
}

void SCF::setOrbitalPrec(double init, double final) {
    this->orbPrec[0] = init;
    this->orbPrec[1] = init;
    this->orbPrec[2] = final;
}

void SCF::adjustPrecision(double error) {
    if (this->orbPrec[0] > 0.0 ) this->orbPrec[0] *= 0.75;
    this->orbPrec[0] = min(10.0*error*error, this->orbPrec[0]);
    this->orbPrec[0] = max(this->orbPrec[0], this->orbPrec[2]);

    this->add.setPrecision(this->orbPrec[0]/10.0);

    TelePrompter::printSeparator(0, '=');
    TelePrompter::printDouble(0, "Current precision", this->orbPrec[0]);
    TelePrompter::printSeparator(0, '-');
    TelePrompter::printDouble(0, "Orbital threshold", this->orbThrs);
    TelePrompter::printDouble(0, "Property threshold", this->propThrs);
    TelePrompter::printSeparator(0, '=', 2);
}

void SCF::resetPrecision() {
    this->orbPrec[0] = this->orbPrec[1];
}

bool SCF::checkConvergence(double err_o, double err_p) const {
    double thrs_o = getOrbitalThreshold();
    double thrs_p = getPropertyThreshold();

    bool conv_o = false;
    bool conv_p = false;
    if (err_o < thrs_o or thrs_o < 0.0) conv_o = true;
    if (err_p < thrs_p or thrs_p < 0.0) conv_p = true;
    return (conv_o and conv_p);
}

bool SCF::needLocalization(int nIter) const {
    bool loc = false;
    if (this->canonical) {
        loc = false;
    } else if (nIter <= 2) {
        loc = true;
    } else if (this->rotation == 0) {
        loc = false;
    } else if (nIter%this->rotation == 0) {
        loc = true;
    }
    return loc;
}

bool SCF::needDiagonalization(int nIter) const {
    bool diag = false;
    if (not this->canonical) {
        diag = false;
    } else if (nIter <= 2) {
        diag = true;
    } else if (this->rotation == 0) {
        diag = false;
    } else if (nIter%this->rotation == 0) {
        diag = true;
    }
    return diag;
}

void SCF::printUpdate(const string &name, double P, double dP) const {
    int oldPrec = TelePrompter::setPrecision(15);
    double p = 1.0;
    if (fabs(P) > MachineZero) {
        p = P;
    }
    double thrs = getPropertyThreshold();
    bool done = (fabs(dP/p) < thrs) or thrs < 0.0;
    printout(0, name);
    printout(0, setw(24) << P);
    TelePrompter::setPrecision(5);
    printout(0, setw(16) << dP);
    println(0, setw(5) << done);
    TelePrompter::setPrecision(oldPrec);
}

double SCF::getUpdate(const vector<double> &vec, int i, bool absPrec) const {
    if (i < 1 or i > vec.size()) MSG_ERROR("Invalid argument");
    double E_i = vec[i-1];
    double E_im1 = 0.0;
    if (i > 1) {
        E_im1 = vec[i-2];
    }
    double E_diff = E_i - E_im1;
    if (not absPrec and fabs(E_i) > MachineZero) {
        E_diff *= 1.0/E_i;
    }
    return E_diff;
}

void SCF::printOrbitals(const VectorXd &epsilon, const OrbitalVector &phi, int flag) const {
    if (mpiOrbRank == 0) {
        TelePrompter::printHeader(0, "Orbitals");
        if (flag == 0) println(0, " Orb    F(i,i)        Error         nNodes  Spin  Occ  Done ");
        if (flag == 1) println(0, " Orb    Norm          Error         nNodes  Spin  Occ  Done ");
        TelePrompter::printSeparator(0, '-');
        int oldprec = TelePrompter::setPrecision(5);
        for (int i = 0; i < phi.size(); i++) {
            const Orbital &phi_i = phi.getOrbital(i);
            printout(0, setw(3) << i);
            printout(0, " " << setw(13) << epsilon(i));
            printout(0, " " << setw(13) << phi_i.getError());
            printout(0, " " << setw(10) << phi_i.getNNodes());
            printout(0, setw(5) << phi_i.printSpin());
            printout(0, setw(5) << phi_i.getOccupancy());
            printout(0, setw(5) << phi_i.isConverged(getOrbitalThreshold()) << endl);
        }

        TelePrompter::printSeparator(0, '-');
        printout(0, " Total error:                    ");
        printout(0, setw(19) << phi.calcTotalError() << "  ");
        printout(0, setw(3) << phi.isConverged(getOrbitalThreshold()) << endl);
        TelePrompter::printSeparator(0, '=', 2);
        TelePrompter::setPrecision(oldprec);
    }
}

void SCF::printConvergence(bool converged) const {
    int iter = this->orbError.size();
    int oldPrec = TelePrompter::getPrecision();
    TelePrompter::printHeader(0, "Convergence rate");
    println(0,"Iter    OrbError       Property                   Update  ");
    TelePrompter::printSeparator(0, '-');
    for (int i = 0; i < iter; i++) {
        double prop_i = this->property[i];
        double propDiff = getUpdate(this->property, i+1, true);
        printout(0, setw(3) << i+1);
        TelePrompter::setPrecision(5);
        printout(0, setw(15) << this->orbError[i]);
        TelePrompter::setPrecision(15);
        printout(0, setw(26) << prop_i);
        TelePrompter::setPrecision(5);
        printout(0, setw(15) << propDiff);
        printout(0, endl);
    }
    TelePrompter::setPrecision(oldPrec);
    TelePrompter::printSeparator(0, '-');
    if (converged) {
        println(0,"                      SCF converged!!!                      ");
    } else {
        println(0,"                   SCF did NOT converge!!!                  ");
    }
    TelePrompter::printSeparator(0, '=', 2);
}

void SCF::printCycle(int nIter) const {
    printout(0, endl << endl);
    printout(0, "#######################");
    printout(0, " SCF cycle " << setw(2) << nIter << " ");
    printout(0, "#######################");
    printout(0, endl << endl << endl);
}

void SCF::printTimer(double t) const {
    int oldPrec = TelePrompter::setPrecision(5);
    printout(0, endl << endl);
    printout(0, "################");
    printout(0, " Wall time: " << t << " sec ");
    printout(0, "################");
    printout(0, endl << endl << endl);
    TelePrompter::setPrecision(oldPrec);
}

void SCF::printMatrix(int level, const MatrixXd &M, const char &name, int pr) const {
    int oldPrec = TelePrompter::setPrecision(pr);
    printout(level, endl);
    printout(level, "----------------------------- ");
    printout(level, name);
    printout(level, " ----------------------------");
    printout(level, endl);
    printout(level, M);
    printout(level, endl);
    printout(level, "------------------------------");
    printout(level, "------------------------------");
    printout(level, endl);
    printout(level, endl);
    TelePrompter::setPrecision(oldPrec);
}
