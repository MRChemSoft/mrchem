#include "HelmholtzOperatorSet.h"
#include "OrbitalVector.h"
#include "OperatorTree.h"
#include "Timer.h"
#include "eigen_disable_warnings.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;
using namespace Eigen;

HelmholtzOperatorSet::HelmholtzOperatorSet(double build, double thrs)
    : threshold(thrs),
      build_prec(build),
      apply(-1.0, MRA->getMaxScale(), true) {
}

void HelmholtzOperatorSet::setup(double prec, const VectorXd &energies) {
    TelePrompter::printHeader(0, "Setting up Helmholtz operators");
    this->apply.setPrecision(prec);

    Timer timer;
    this->operIdx.clear();
    this->lambda.clear();
    if (mpiOrbSize > 1) this->clear();
    for (int i = 0; i < energies.size(); i++) {
        double energy = energies(i);
        if (energy > 0.0 ) energy = -0.5;
        int idx = initHelmholtzOperator(energy, i);
	
        this->operIdx.push_back(idx);
	if (mpiOrbSize == 1) {
	    double mu = getOperator(i).getMu();
	    this->lambda.push_back(-0.5*mu*mu);
	} else { 
	    this->lambda.push_back(energy);
        }
    }
    if (mpiOrbSize == 1) clearUnused();
    timer.stop();
    TelePrompter::printFooter(0, timer, 2);
}

int HelmholtzOperatorSet::initHelmholtzOperator(double energy, int i) {
    if (energy > 0.0) MSG_ERROR("Complex Helmholtz not available");
    double mu = sqrt(-2.0*energy);
    if (mpiOrbSize == 1) {
	for (int j = 0; j < this->operators.size(); j++) {
	    double mu_j = this->operators[j]->getMu();
	    double muDiff = mu - mu_j;
	    double relDiff = muDiff/mu;
	    if (fabs(relDiff) < this->threshold and false) {
		double l = -0.5*mu_j*mu_j;
		TelePrompter::printDouble(0, "Re-using operator with lambda", l);
		return j;
	    }
	}
    }
    HelmholtzOperator *oper = 0;
    if( i%mpiOrbSize == mpiOrbRank) {
	TelePrompter::printDouble(0, "Creating operator with lambda", energy);	
	oper = new HelmholtzOperator(*MRA, mu, this->build_prec);
    }
    this->operators.push_back(oper);

    return this->operators.size() - 1;
}

void HelmholtzOperatorSet::clear() {
    int nOper = this->operators.size();
    for (int j = 0; j < nOper; j++) {
        if (this->operators[j] != 0) {
            delete this->operators[j];
            this->operators[j] = 0;
        }
    }
    this->operators.clear();
    this->operIdx.clear();
    this->lambda.clear();
}

void HelmholtzOperatorSet::clearUnused() {
    vector<HelmholtzOperator *> tmp;
    int nIdx = this->operIdx.size();
    int nOper = this->operators.size();
    VectorXd nUsed = VectorXd::Zero(nOper);
    for (int i = 0; i < nIdx; i++) {
        int n = this->operIdx[i];
        nUsed(n) = nUsed(n) + 1;
    }
    for (int i = 0; i < nOper; i++) {
        if (nUsed(i) == 0) {
            delete this->operators[i];
        } else {
            tmp.push_back(this->operators[i]);
            int newIdx = tmp.size() - 1;
            for (int j = 0; j < nIdx; j++) {
                if (operIdx[j] == i) operIdx[j] = newIdx;
            }
        }
        this->operators[i] = 0;
    }
    this->operators.clear();
    for (int i = 0; i < tmp.size(); i++) {
        this->operators.push_back(tmp[i]);
        tmp[i] = 0;
    }
}

HelmholtzOperator& HelmholtzOperatorSet::getOperator(int i) {
    int idx = this->operIdx[i];
    if (idx < 0 or idx >= operators.size()) MSG_ERROR("Invalid operator index");
    HelmholtzOperator *oper = this->operators[idx];
    if (oper == 0) MSG_ERROR("Operator null pointer");
    return *oper;
}

VectorXd HelmholtzOperatorSet::getLambda() const {
    int nLambda = this->lambda.size();
    VectorXd L = VectorXd::Zero(nLambda);
    for (int i = 0; i < nLambda; i++) {
        L(i) = this->lambda[i];
    }
    return L;
}

/** Prints the number of trees and nodes kept in the operator set */

int HelmholtzOperatorSet::printTreeSizes() const {
    int totNodes = 0;
    int totTrees = 0;
    int nOperators = this->operators.size();
    int k = MRA->getOrder();
    int kp1 = k + 1;
    for (int i = 0; i < nOperators; i++) {
	if(this->operators[i] != 0){
	    int nTrees = this->operators[i]->size();
	    for (int j = 0; j < nTrees; j++) {
		totNodes += this->operators[i]->getComponent(j).getNNodes();
	    }
	    totTrees += nTrees;
	}
    }
    println(0, " HelmholtzOperators " << totTrees <<" trees, " << totNodes<<" nodes, "<<nOperators<<" operators, total"<< totNodes*49*4*8/1024/1024<<" MB");

    return totNodes;
}


void HelmholtzOperatorSet::operator()(int i, Orbital &out, Orbital &inp) {
    if (out.hasReal()) MSG_ERROR("Orbital not empty");
    if (out.hasImag()) MSG_ERROR("Orbital not empty");

    HelmholtzOperator &H_i = getOperator(i);

    if (inp.hasReal()) {
        out.allocReal();
        this->apply(out.real(), H_i, inp.real());
    }
    if (inp.hasImag()) {
        out.allocImag();
        this->apply(out.imag(), H_i, inp.imag());
    }
}

void HelmholtzOperatorSet::operator()(OrbitalVector &out, OrbitalVector &inp) {
    TelePrompter::printHeader(0, "Applying Helmholtz operators");
    println(0, " Orb  RealNorm     ImagNorm      nNodes     Error   Timing   ");
    TelePrompter::printSeparator(0, '-');
    int oldprec = TelePrompter::setPrecision(5);

    HelmholtzOperatorSet &H = *this;

    Timer tottimer;
    for (int i = mpiOrbRank; i < inp.size(); i += mpiOrbSize) {
        Timer timer;
        Orbital &out_i = out.getOrbital(i);
        Orbital &inp_i = inp.getOrbital(i);
        H(i, out_i, inp_i);

        int rNodes = out_i.getNNodes(QMFunction::Real);
        int iNodes = out_i.getNNodes(QMFunction::Imag);
        double norm = sqrt(out_i.getSquareNorm());
        double dNorm = fabs(norm-1.0);
        double real_norm = sqrt(out_i.getSquareNorm(QMFunction::Real));
        double imag_norm = sqrt(out_i.getSquareNorm(QMFunction::Imag));

        timer.stop();
        TelePrompter::setPrecision(5);
        cout << setw(3) << i;
        cout << " " << setw(12) << real_norm;
        cout << " " << setw(12) << imag_norm;
        TelePrompter::setPrecision(1);
        cout << " " << setw(5) << rNodes;
        cout << " " << setw(5) << iNodes;
        cout << " " << setw(8) << dNorm;
        cout << setw(9) << timer.getWallTime() << endl;
    }

#ifdef HAVE_MPI
    // barrier to align printing
    MPI_Barrier(mpiCommOrb);
#endif

    tottimer.stop();
    TelePrompter::printFooter(0, tottimer, 2);
    TelePrompter::setPrecision(oldprec);
}
