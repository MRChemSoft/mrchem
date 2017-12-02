#include "FockOperator.h"
#include "KineticOperator.h"
#include "NuclearPotential.h"
#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "XCOperator.h"
#include "OrbitalAdder.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "SCFEnergy.h"
#include "MathUtils.h"
#include "Timer.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA
extern OrbitalVector workOrbVec;

using namespace std;
using namespace Eigen;

FockOperator::FockOperator(KineticOperator *t,
                           NuclearPotential *v,
                           CoulombOperator *j,
                           ExchangeOperator *k,
                           XCOperator *xc)
    : QMOperator(MRA->getMaxScale()),
      T(t),
      V(v),
      J(j),
      K(k),
      XC(xc),
      H_1(0) {
}

FockOperator::~FockOperator() {
    this->T = 0;
    this->V = 0;
    this->J = 0;
    this->K = 0;
    this->XC = 0;
    this->H_1 = 0;
}

void FockOperator::setup(double prec) {
    Timer timer;
    TelePrompter::printHeader(0, "Setting up Fock operator");
    TelePrompter::printDouble(0, "Precision", prec);
    TelePrompter::printSeparator(0, '-');
    if (this->T != 0) this->T->setup(prec);
    if (this->V != 0) this->V->setup(prec);
    if (this->J != 0) this->J->setup(prec);
    if (this->K != 0) this->K->setup(prec);
    if (this->XC != 0) this->XC->setup(prec);
    if (this->H_1 != 0) this->H_1->setup(prec);
    timer.stop();
    TelePrompter::printFooter(0, timer, 2);
}

void FockOperator::clear() {
    if (this->T != 0) this->T->clear();
    if (this->V != 0) this->V->clear();
    if (this->J != 0) this->J->clear();
    if (this->K != 0) this->K->clear();
    if (this->XC != 0) this->XC->clear();
    if (this->H_1 != 0) this->H_1->clear();
}

void FockOperator::rotate(MatrixXd &U) {
    if (this->K != 0) this->K->rotate(U);
}

Orbital* FockOperator::operator() (Orbital &orb_p) {
    OrbitalAdder add(1.0, this->max_scale);

    vector<Orbital *> orbs;
    if (this->T != 0) orbs.push_back((*this->T)(orb_p));
    if (this->V != 0) orbs.push_back((*this->V)(orb_p));
    if (this->J != 0) orbs.push_back((*this->J)(orb_p));
    if (this->K != 0) orbs.push_back((*this->K)(orb_p));
    if (this->XC != 0) orbs.push_back((*this->XC)(orb_p));
    if (this->H_1 != 0) orbs.push_back((*this->H_1)(orb_p));

    vector<complex<double> > coefs;
    for (int i = 0; i < orbs.size(); i++) coefs.push_back(1.0);

    Orbital *result = new Orbital(orb_p);
    add(*result, coefs, orbs, true);

    for (int n = 0; n < orbs.size(); n++) {
        if (orbs[n] != 0) delete orbs[n];
        orbs[n] = 0;
    }
    return result;
}

Orbital* FockOperator::adjoint(Orbital &orb_p) {
    NOT_IMPLEMENTED_ABORT;
}

double FockOperator::operator() (Orbital &orb_i, Orbital &orb_j) {
    double result = 0.0;
    if (this->T != 0) result += (*this->T)(orb_i, orb_j);
    if (this->V != 0) result += (*this->V)(orb_i, orb_j);
    if (this->J != 0) result += (*this->J)(orb_i, orb_j);
    if (this->K != 0) result += (*this->K)(orb_i, orb_j);
    if (this->XC != 0) result += (*this->XC)(orb_i, orb_j);
    if (this->H_1 != 0) result += (*this->H_1)(orb_i, orb_j);
    return result;
}

double FockOperator::adjoint(Orbital &orb_i, Orbital &orb_j) {
    NOT_IMPLEMENTED_ABORT;
}


MatrixXd FockOperator::operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    MatrixXd result = MatrixXd::Zero(Ni,Nj);
    Timer tot_t;
    TelePrompter::printHeader(0, "Calculating Fock matrix");

#ifdef HAVE_MPI
    OrbitalVector orbVecChunk_i(0); //to store adresses of own i_orbs
    OrbitalVector orbVecChunk_j(0); //to store adresses of own j_orbs
    OrbitalVector rcvOrbs(0);       //to store adresses of received orbitals
    vector<int> orbsIx;             //to store own orbital indices
    int rcvOrbsIx[workOrbVecSize];  //to store received orbital indices

    //make vector with adresses of own orbitals
    for (int Ix = mpiOrbRank; Ix < Ni; Ix += mpiOrbSize) {
	orbVecChunk_i.push_back(i_orbs.getOrbital(Ix));//i orbitals
	orbsIx.push_back(Ix);
    }
    for (int Ix = mpiOrbRank; Ix < Nj; Ix += mpiOrbSize) orbVecChunk_j.push_back(j_orbs.getOrbital(Ix));//j orbitals
    //need to pad orbVecChunk_j so that all have same size
    if (orbVecChunk_j.size() < (Nj+mpiOrbSize-1)/mpiOrbSize) orbVecChunk_j.push_back(j_orbs.getOrbital(mpiOrbRank));

    for (int iter = 0; iter >= 0; iter++) {
	//get a new chunk from other processes
	//NB: should not use directly workorbvec as rcvOrbs, because they may 
	//contain own orbitals, and these can be overwritten
	orbVecChunk_i.getOrbVecChunk(orbsIx, rcvOrbs, rcvOrbsIx, Ni, iter);

	//Only one process does the computations. j orbitals always local
	MatrixXd resultChunk = MatrixXd::Zero(rcvOrbs.size(),orbVecChunk_j.size());
      
	if (this->T != 0) resultChunk = (*this->T)(rcvOrbs,orbVecChunk_j);
	if (this->V != 0) resultChunk += (*this->V)(rcvOrbs,orbVecChunk_j);
	if (this->J != 0) resultChunk += (*this->J)(rcvOrbs,orbVecChunk_j);
	if (this->K != 0) {
	    resultChunk += (*this->K)(rcvOrbs,orbVecChunk_j);
	    if (rcvOrbs.size() == 0 and iter>=0) {
		//we must still go through operator to send own orbitals to others. Just make a fake operation!
		rcvOrbs.push_back(i_orbs.getOrbital(mpiOrbRank));//i orbitals
		MatrixXd resultChunk_notused = MatrixXd::Zero(rcvOrbs.size(),orbVecChunk_j.size());
		resultChunk_notused = (*this->K)(rcvOrbs,orbVecChunk_j);
	    }
	}
	if (this->XC != 0) resultChunk += (*this->XC)(rcvOrbs,orbVecChunk_j);

	//copy results into final matrix
	int j = 0;
	for (int Jx = mpiOrbRank; Jx < Nj; Jx += mpiOrbSize) {
	    for (int ix = 0; ix < rcvOrbs.size(); ix++) {
		result(rcvOrbsIx[ix],Jx) += resultChunk(ix,j);
	    }
	    j++;
	}
	rcvOrbs.clearVec(false);//reset to zero size orbital vector
    }

    //clear orbital vector adresses. NB: only references and metadata must be deleted, not the trees in orbitals
    orbVecChunk_i.clearVec(false);
    orbVecChunk_j.clearVec(false);
    workOrbVec.clear();

    //combine results from all processes
    MPI_Allreduce(MPI_IN_PLACE, &result(0,0), Ni*Nj,
                  MPI_DOUBLE, MPI_SUM, mpiCommOrb);

#else
    if (this->T != 0) {
        Timer timer;
        result += (*this->T)(i_orbs, j_orbs);
        timer.stop();
        TelePrompter::printDouble(0, "Kinetic matrix", timer.getWallTime());
    }
    if (this->V != 0) {
        Timer timer;
        result += (*this->V)(i_orbs, j_orbs);
        timer.stop();
        TelePrompter::printDouble(0, "Nuclear potential matrix", timer.getWallTime());
    }
    if (this->J != 0) {
        Timer timer;
        result += (*this->J)(i_orbs, j_orbs);
        timer.stop();
        TelePrompter::printDouble(0, "Coulomb matrix", timer.getWallTime());
    }
    if (this->K != 0) {
        Timer timer;
        result += (*this->K)(i_orbs, j_orbs);
        timer.stop();
        TelePrompter::printDouble(0, "Exchange matrix", timer.getWallTime());
    }
    if (this->XC != 0) {
        Timer timer;
        result += (*this->XC)(i_orbs, j_orbs);
        timer.stop();
        TelePrompter::printDouble(0, "Exchange-Correlation matrix", timer.getWallTime());
    }
    if (this->H_1 != 0) {
        Timer timer;
        result += (*this->H_1)(i_orbs, j_orbs);
        timer.stop();
        TelePrompter::printDouble(0, "Perturbation matrix", timer.getWallTime());
    }
#endif
    tot_t.stop();
    TelePrompter::printFooter(0, tot_t, 2);

    return result;
}

MatrixXd FockOperator::adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    if(mpiOrbSize>1)cout<<"ERROR"<<endl;
    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    MatrixXd result = MatrixXd::Zero(Ni,Nj);
    if (this->T != 0) result += (*this->T).adjoint(i_orbs, j_orbs);
    if (this->V != 0) result += (*this->V).adjoint(i_orbs, j_orbs);
    if (this->J != 0) result += (*this->J).adjoint(i_orbs, j_orbs);
    if (this->K != 0) result += (*this->K).adjoint(i_orbs, j_orbs);
    if (this->XC != 0) result += (*this->XC).adjoint(i_orbs, j_orbs);
    if (this->H_1 != 0) result += (*this->H_1).adjoint(i_orbs, j_orbs);
    return result;
}

Orbital* FockOperator::applyPotential(Orbital &orb_p) {
    OrbitalAdder add(1.0, this->max_scale);

    vector<Orbital *> orbs;
    if (this->V != 0) orbs.push_back((*this->V)(orb_p));
    if (this->J != 0) orbs.push_back((*this->J)(orb_p));
    if (this->K != 0) orbs.push_back((*this->K)(orb_p));
    if (this->XC != 0) orbs.push_back((*this->XC)(orb_p));

    vector<complex<double> > coefs;
    for (int i = 0; i < orbs.size(); i++) coefs.push_back(1.0);

    Timer timer;
    Orbital *result = new Orbital(orb_p);
    add(*result, coefs, orbs, true);
    timer.stop();
    double time = timer.getWallTime();
    int nNodes = result->getNNodes();
    if(mpiOrbSize==1)TelePrompter::printTree(1, "Sum potential operator", nNodes, time);

    for (int n = 0; n < orbs.size(); n++) {
        delete orbs[n];
        orbs[n] = 0;
    }
    return result;
}

double FockOperator::applyPotential(Orbital &orb_i, Orbital &orb_j) {
    double result = 0.0;
    if (this->V != 0) result += (*this->V)(orb_i, orb_j);
    if (this->J != 0) result += (*this->J)(orb_i, orb_j);
    if (this->K != 0) result += (*this->K)(orb_i, orb_j);
    if (this->XC != 0) result += (*this->XC)(orb_i, orb_j);
    return result;
}

MatrixXd FockOperator::applyPotential(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    int nOrbs = i_orbs.size();
    MatrixXd result = MatrixXd::Zero(nOrbs,nOrbs);
    if (this->V != 0) result += (*this->V)(i_orbs, j_orbs);
    if (this->J != 0) result += (*this->J)(i_orbs, j_orbs);
    if (this->K != 0) result += (*this->K)(i_orbs, j_orbs);
    if (this->XC != 0) result += (*this->XC)(i_orbs, j_orbs);
    return result;
}

Orbital* FockOperator::applyAdjointPotential(Orbital &orb_p) {
    NOT_IMPLEMENTED_ABORT;
}

double FockOperator::applyAdjointPotential(Orbital &orb_i, Orbital &orb_j) {
    NOT_IMPLEMENTED_ABORT;
}

MatrixXd FockOperator::applyAdjointPotential(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    NOT_IMPLEMENTED_ABORT;
}

SCFEnergy FockOperator::trace(OrbitalVector &phi, MatrixXd &F) {
    double E_nuc = 0.0;
    double E_el = 0.0;
    double E_orb = 0.0;
    double E_kin = 0.0;
    double E_en = 0.0;
    double E_ee = 0.0;
    double E_xc = 0.0;
    double E_x = 0.0;

    // Nuclear part
    if (this->V != 0) {
        Nuclei &nucs = this->V->getNuclei();
        int nNucs = nucs.size();
        for (int i = 0; i < nNucs; i++) {
            for (int j = i+1; j < nNucs; j++) {
                const Nucleus &nuc_i = nucs[i];
                const Nucleus &nuc_j = nucs[j];
                double Z_i = nuc_i.getCharge();
                double Z_j = nuc_j.getCharge();
                const double *R_i = nuc_i.getCoord();
                const double *R_j = nuc_j.getCoord();
                double r_ij = MathUtils::calcDistance(3, R_i, R_j);
                E_nuc += (Z_i*Z_j)/r_ij;
            }
        }
    }

    // Electronic part
    double E_xc2 = 0.0;
    if (this->XC != 0) E_xc = this->XC->getEnergy();
    for (int i = 0; i < phi.size(); i++) {
        if (i%mpiOrbSize == mpiOrbRank) {
            Orbital &phi_i = phi.getOrbital(i);
            double occ = (double) phi_i.getOccupancy();
            double e_i = occ*F(i,i);
            E_orb += e_i;

            if (this->V != 0) {
                println(2, "\n<" << i << "|V_nuc|" << i << ">");
                E_en += occ*(*this->V)(phi_i,phi_i);
            }
            if (this->J != 0) {
                println(2, "\n<" << i << "|J|" << i << ">");
                E_ee += 0.5*occ*(*this->J)(phi_i,phi_i);
            }
            if (this->K != 0) {
                println(2, "\n<" << i << "|K|" << i << ">");
                E_x += 0.5*occ*(*this->K)(phi_i,phi_i);
            }
            if (this->XC != 0) {
                println(2, "\n<" << i << "|V_xc|" << i << ">");
                E_xc2 += occ*(*this->XC)(phi_i,phi_i);
            }
        }
    }

#ifdef HAVE_MPI
    double tmp[5];
    tmp[0] = E_orb;
    tmp[1] = E_en;
    tmp[2] = E_ee;
    tmp[3] = E_x;
    tmp[4] = E_xc2;
    MPI_Allreduce(MPI_IN_PLACE,tmp, 5, MPI_DOUBLE, MPI_SUM, mpiCommOrb);
    E_orb = tmp[0];
    E_en  = tmp[1];
    E_ee  = tmp[2];
    E_x   = tmp[3];
    E_xc2 = tmp[4];
#endif

    double E_eex = E_ee + E_x;
    double E_orbxc2 = E_orb - E_xc2;
    E_kin = E_orbxc2 - 2.0*E_eex - E_en;
    E_el = E_orbxc2 - E_eex + E_xc;

    return SCFEnergy(E_nuc, E_el, E_orb, E_kin, E_en, E_ee, E_xc, E_x);
}
