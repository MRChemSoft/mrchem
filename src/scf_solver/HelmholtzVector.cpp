#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"

#include "HelmholtzVector.h"
#include "Orbital.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

HelmholtzVector::HelmholtzVector(double build, double thrs)
        : threshold(thrs),
          build_prec(build),
          apply_prec(-1.0) {
}

void HelmholtzVector::setup(double prec, const DoubleVector &energies) {
    Printer::printHeader(0, "Setting up Helmholtz operators");
    this->apply_prec = prec;

    Timer timer;
    this->oper_idx.clear();
    this->lambda.clear();
    if (mpi::orb_size > 1) this->clear();
    for (int i = 0; i < energies.size(); i++) {
        double energy = energies(i);
        if (energy > 0.0 ) energy = -0.5;
        int idx = initHelmholtzOperator(energy, i);
        this->oper_idx.push_back(idx);
        if (mpi::orb_size == 1) {
            double mu = (*this)[i].getMu();
            this->lambda.push_back(-0.5*mu*mu);
        } else {
            this->lambda.push_back(energy);
        }
    }
    if (mpi::orb_size == 1) clearUnused();

    timer.stop();
    Printer::printFooter(0, timer, 2);
}

int HelmholtzVector::initHelmholtzOperator(double energy, int i) {
    if (energy > 0.0) MSG_ERROR("Complex Helmholtz not available");
    double mu = std::sqrt(-2.0*energy);
    if (mpi::orb_size == 1) {
        for (int j = 0; j < this->operators.size(); j++) {
            double mu_j = this->operators[j]->getMu();
            double muDiff = mu - mu_j;
            double relDiff = muDiff/mu;
            if (std::abs(relDiff) < this->threshold and false) {
                double l = -0.5*mu_j*mu_j;
                Printer::printDouble(0, "Re-using operator with lambda", l, 5);
                return j;
            }
        }
    }
    mrcpp::HelmholtzOperator *oper = 0;
    if( i%mpi::orb_size == mpi::orb_rank) {
        Printer::printDouble(0, "Creating operator with lambda", energy, 5);
        oper = new mrcpp::HelmholtzOperator(*MRA, mu, this->build_prec);
    }
    this->operators.push_back(oper);

    return this->operators.size() - 1;
}

void HelmholtzVector::clear() {
    int nOper = this->operators.size();
    for (int j = 0; j < nOper; j++) {
        if (this->operators[j] != 0) {
            delete this->operators[j];
            this->operators[j] = 0;
        }
    }
    this->operators.clear();
    this->oper_idx.clear();
    this->lambda.clear();
}

void HelmholtzVector::clearUnused() {
    std::vector<mrcpp::HelmholtzOperator *> tmp;
    int nIdx = this->oper_idx.size();
    int nOper = this->operators.size();
    IntVector nUsed = IntVector::Zero(nOper);
    for (int i = 0; i < nIdx; i++) {
        int n = this->oper_idx[i];
        nUsed(n) = nUsed(n) + 1;
    }
    for (int i = 0; i < nOper; i++) {
        if (nUsed(i) == 0) {
            delete this->operators[i];
        } else {
            tmp.push_back(this->operators[i]);
            int newIdx = tmp.size() - 1;
            for (int j = 0; j < nIdx; j++) {
                if (this->oper_idx[j] == i) this->oper_idx[j] = newIdx;
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

mrcpp::HelmholtzOperator& HelmholtzVector::operator[](int i) {
    int idx = this->oper_idx[i];
    if (idx < 0 or idx >= operators.size()) MSG_ERROR("Invalid operator index");
    mrcpp::HelmholtzOperator *oper = this->operators[idx];
    if (oper == 0) MSG_ERROR("Operator null pointer");
    return *oper;
}

const mrcpp::HelmholtzOperator& HelmholtzVector::operator[](int i) const {
    int idx = this->oper_idx[i];
    if (idx < 0 or idx >= operators.size()) MSG_ERROR("Invalid operator index");
    const mrcpp::HelmholtzOperator *oper = this->operators[idx];
    if (oper == 0) MSG_ERROR("Operator null pointer");
    return *oper;
}

DoubleVector HelmholtzVector::getLambda() const {
    int nLambda = this->lambda.size();
    DoubleVector L = DoubleVector::Zero(nLambda);
    for (int i = 0; i < nLambda; i++) {
        L(i) = this->lambda[i];
    }
    return L;
}

/** Prints the number of trees and nodes kept in the operator set */

int HelmholtzVector::printTreeSizes() const {
    int totNodes = 0;
    int totTrees = 0;
    int nOperators = this->operators.size();
    for (int i = 0; i < nOperators; i++) {
        if(this->operators[i] != 0){
            int nTrees = this->operators[i]->size();
            for (int j = 0; j < nTrees; j++) {
                totNodes += this->operators[i]->getComponent(j).getNNodes();
            }
            totTrees += nTrees;
        }
    }
    return totNodes;
}


Orbital HelmholtzVector::operator()(int i, Orbital inp) {
    mrcpp::HelmholtzOperator &H_i = (*this)[i];

    Orbital out = inp.paramCopy();
    if (inp.hasReal()) {
        out.alloc(NUMBER::Real);
        mrcpp::apply(this->apply_prec, out.real(), H_i, inp.real());
    }
    if (inp.hasImag()) {
        out.alloc(NUMBER::Imag);
        mrcpp::apply(this->apply_prec, out.imag(), H_i, inp.imag());
        if (inp.conjugate()) out.imag().rescale(-1.0);
    }
    return out;
}

OrbitalVector HelmholtzVector::operator()(OrbitalVector &inp) {
    Printer::printHeader(0, "Applying Helmholtz operators");
    println(0, " Orb  RealNorm     ImagNorm      nNodes     Error   Timing   ");
    Printer::printSeparator(0, '-');
    int oldprec = Printer::setPrecision(5);

    OrbitalVector out = orbital::param_copy(inp);
    HelmholtzVector &H = (*this);

    Timer tottimer;
    for (int i = 0; i < inp.size(); i++) {
        if (not mpi::my_orb(inp[i])) continue;

        Timer timer;
        out[i] = H(i, inp[i]);

        int rNodes = out[i].getNNodes(NUMBER::Real);
        int iNodes = out[i].getNNodes(NUMBER::Imag);
        double dNorm = std::abs(out[i].norm() - 1.0);
        double rNorm = 0.0;
        double iNorm = 0.0;
        if (out[i].hasReal()) rNorm = std::sqrt(out[i].real().getSquareNorm());
        if (out[i].hasImag()) iNorm = std::sqrt(out[i].imag().getSquareNorm());

        timer.stop();
        Printer::setPrecision(5);
        printout(0, std::setw(3) << i);
        printout(0, " " << std::setw(12) << rNorm);
        printout(0, " " << std::setw(12) << iNorm);
        Printer::setPrecision(1);
        printout(0, " " << std::setw(5) << rNodes);
        printout(0, " " << std::setw(5) << iNodes);
        printout(0, " " << std::setw(8) << dNorm);
        printout(0, std::setw(9) << timer.getWallTime() << std::endl);
    }

#ifdef HAVE_MPI
    MPI_Barrier(mpi::comm_orb); // barrier to align printing
#endif

    tottimer.stop();
    Printer::printFooter(0, tottimer, 2);
    Printer::setPrecision(oldprec);
    return out;
}

} //namespace mrchem
