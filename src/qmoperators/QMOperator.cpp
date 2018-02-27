#include "QMOperator.h"
#include "Orbital.h"
#include "OrbitalVector.h"
#include "Timer.h"

using namespace Eigen;
using namespace std;

void QMOperator::setApplyPrec(double prec) {
    if (this->apply_prec < 0.0) { 
        this->apply_prec = prec;
    } else if (not IS_EQUAL(prec, this->apply_prec)) {
        MSG_ERROR("Clear operator before setup with different prec!");
    }
}

double QMOperator::operator() (Orbital &phi_i, Orbital &phi_j) {
    QMOperator &O = *this;
    Orbital *Ophi_j = O(phi_j);
    complex<double> result = phi_i.dot(*Ophi_j);
    delete Ophi_j;
    if (result.imag() > MachineZero) MSG_ERROR("Should be real");
    return result.real();
}

double QMOperator::adjoint(Orbital &phi_i, Orbital &phi_j) {
    NOT_IMPLEMENTED_ABORT;
}

MatrixXd QMOperator::operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    if (mpiOrbSize == 1) TelePrompter::printHeader(1, "Compute Matrix Element");
    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    MatrixXcd result = MatrixXcd::Zero(Ni, Nj);

    Timer timer;
    QMOperator &O = *this;
    for (int j = 0; j < Nj; j++) {
        Orbital &phi_j = j_orbs.getOrbital(j);
        Orbital *Ophi_j = O(phi_j);
        for (int i = 0; i <  Ni; i++) {
            Orbital &phi_i = i_orbs.getOrbital(i);
            result(i,j) = phi_i.dot(*Ophi_j);
        }
        delete Ophi_j;
    }
    if (mpiOrbSize == 1) timer.stop();
    if (mpiOrbSize == 1) TelePrompter::printFooter(1, timer, 2);

    if (result.imag().norm() > MachineZero) MSG_ERROR("Should be real");
    return result.real();
}

MatrixXd QMOperator::adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    NOT_IMPLEMENTED_ABORT;
}

