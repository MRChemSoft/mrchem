#include "QMOperatorExp.h"
#include "Orbital.h"
#include "OrbitalVector.h"
#include "OrbitalAdder.h"
#include "TelePrompter.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace Eigen;
using namespace std;

QMOperatorExp& QMOperatorExp::operator+=(QMOperator &O) {
    this->coefs.push_back(1.0);
    vector<QMOperator *> tmp;
    tmp.push_back(&O);
    this->oper_exp.push_back(tmp);
    return *this;
}

QMOperatorExp& QMOperatorExp::operator=(QMOperator &O) {
    this->clear();
    this->coefs.push_back(1.0);
    vector<QMOperator *> tmp;
    tmp.push_back(&O);
    this->oper_exp.push_back(tmp);
    return *this;
}

QMOperatorExp& QMOperatorExp::operator+=(const QMOperatorExp &O) {
    if (this != &O) {
        for (int i = 0; i < O.coefs.size(); i++) {
            this->coefs.push_back(O.coefs[i]);
        }
        for (int i = 0; i < O.oper_exp.size(); i++) {
            this->oper_exp.push_back(O.oper_exp[i]);
        }
    }
    return *this;
}

QMOperatorExp& QMOperatorExp::operator=(const QMOperatorExp &O) {
    if (this != &O) {
        this->coefs = O.coefs;
        this->oper_exp = O.oper_exp;
    }
    return *this;
}

void QMOperatorExp::setup(double prec) {
    for (int i = 0; i < this->oper_exp.size(); i++) {
        for (int j = 0; j < this->oper_exp[i].size(); j++) {
            this->oper_exp[i][j]->setup(prec);
        }
    }
}
void QMOperatorExp::clear() {
    for (int i = 0; i < this->oper_exp.size(); i++) {
        for (int j = 0; j < this->oper_exp[i].size(); j++) {
            this->oper_exp[i][j]->clear();
        }
    }
}

Orbital* QMOperatorExp::operator() (Orbital &phi_p) {
    OrbitalAdder add(-1.0, MRA->getMaxScale());

    vector<Orbital *> orbs;
    for (int i = 0; i < this->oper_exp.size(); i++) {
        if (this->oper_exp[i][0] == 0) MSG_FATAL("Invalid oper exp");
        QMOperator &O_i0 = *this->oper_exp[i][0];
        Orbital *tmp_1 = O_i0(phi_p);
        for (int j = 1; j < this->oper_exp[i].size(); j++) {
            if (this->oper_exp[i][j] == 0) MSG_FATAL("Invalid oper exp");
            QMOperator &O_ij = *this->oper_exp[i][j];
            Orbital *tmp_2 = O_ij(*tmp_1);
            delete tmp_1;
            tmp_1 = tmp_2;
        }
        orbs.push_back(tmp_1);
    }
    Orbital *result = new Orbital(phi_p);
    add(*result, this->coefs, orbs, true);

    for (int n = 0; n < orbs.size(); n++) {
        if (orbs[n] != 0) delete orbs[n];
        orbs[n] = 0;
    }
    return result;
}

Orbital* QMOperatorExp::adjoint(Orbital &phi_p) {
    NOT_IMPLEMENTED_ABORT;
}

double QMOperatorExp::operator() (Orbital &phi_i, Orbital &phi_j) {
    QMOperatorExp &O = *this;
    Orbital *Ophi_j = O(phi_j);
    complex<double> result = phi_i.dot(*Ophi_j);
    delete Ophi_j;
    if (result.imag() > MachineZero) MSG_ERROR("Should be real");
    return result.real();
}

double QMOperatorExp::adjoint(Orbital &phi_i, Orbital &phi_j) {
    NOT_IMPLEMENTED_ABORT;
}

MatrixXd QMOperatorExp::operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    QMOperatorExp &O = *this;

    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    MatrixXcd result = MatrixXcd::Zero(Ni, Nj);

    for (int j = 0; j < Nj; j++) {
        Orbital &phi_j = j_orbs.getOrbital(j);
        Orbital *Ophi_j = O(phi_j);
        for (int i = 0; i <  Ni; i++) {
            Orbital &phi_i = i_orbs.getOrbital(i);
            result(i,j) = phi_i.dot(*Ophi_j);
        }
        delete Ophi_j;
    }

    if (result.imag().norm() > MachineZero) MSG_ERROR("Should be real");
    return result.real();
}

MatrixXd QMOperatorExp::adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    NOT_IMPLEMENTED_ABORT;
}

double QMOperatorExp::trace(OrbitalVector &phi, OrbitalVector &x, OrbitalVector &y) {
    QMOperatorExp &O = *this;

    double result = 0.0;
    for (int i = 0; i < phi.size(); i++) {
        if (i%mpiOrbSize == mpiOrbRank) {
	    Orbital &phi_i = phi.getOrbital(i);
	    Orbital &x_i = x.getOrbital(i);
	    Orbital &y_i = y.getOrbital(i);
	    
	    double eta_i = (double) phi_i.getOccupancy();
	    double result_1 = O(phi_i, x_i);
	    double result_2 = O(y_i, phi_i);
	    result += eta_i*(result_1 + result_2);
        }
    }
#ifdef HAVE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, mpiCommOrb);
#endif

    return result;
}

double QMOperatorExp::trace(OrbitalVector &phi) {
    QMOperatorExp &O = *this;

    double result = 0.0;
    for (int i = 0; i < phi.size(); i++) {
	if (i%mpiOrbSize == mpiOrbRank) {
	    Orbital &phi_i = phi.getOrbital(i);
	    double eta_i = (double) phi_i.getOccupancy();
	    result += eta_i*O(phi_i, phi_i);
        }
    }
#ifdef HAVE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, mpiCommOrb);
#endif
    return result;
}
