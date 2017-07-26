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
  if (MPI_rank == 0) {
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

/** Computes a new set of orbitals by application of the Helmholtz operator
 *
 * Requires orbitals in the orbitals pointer as well as corresponding operators
 * in the Helmholtz set and produces new non-orthogonal orbitals which are
 * stored locally in the newOrbs pointer. Does not compute orbital differences,
 * so the orbUpdates set is empty after this routine. Errors are estimated by
 * the norms of the new orbitals (deviation from one).
 */
void SCF::applyHelmholtzOperators(OrbitalVector &phi_np1,
                                  MatrixXd &F_n,
                                  OrbitalVector &phi_n,
                                  bool adjoint) {
    TelePrompter::printHeader(0, "Applying Helmholtz Operators");
    println(0, " Orb  RealNorm     ImagNorm      nNodes     Error   Timing   ");
    TelePrompter::printSeparator(0, '-');
    int oldprec = TelePrompter::setPrecision(5);

    Timer timer;
    HelmholtzOperatorSet &H = *this->helmholtz;
    H.setPrecision(getOrbitalPrecision());
    phi_np1.clear();
    for (int i = 0; i < phi_n.size(); i++) {
        Timer timer;
        Orbital &nPhi_i = phi_n.getOrbital(i);
        Orbital &np1Phi_i = phi_np1.getOrbital(i);

	if (i%MPI_size == MPI_rank) {
	    //in charge for this orbital
	    Orbital *arg_i = getHelmholtzArgument(i, F_n, phi_n, adjoint);
	    H(i, np1Phi_i, *arg_i);
	    delete arg_i;
	  
        int rNodes = np1Phi_i.getNNodes(QMFunction::Real);
        int iNodes = np1Phi_i.getNNodes(QMFunction::Imag);
	    double norm_n = sqrt(nPhi_i.getSquareNorm());
	    double norm_np1 = sqrt(np1Phi_i.getSquareNorm());
	    double dNorm_n = fabs(norm_np1-norm_n);
        double real_norm = sqrt(np1Phi_i.getSquareNorm(QMFunction::Real));
        double imag_norm = sqrt(np1Phi_i.getSquareNorm(QMFunction::Imag));

	    timer.stop();
            TelePrompter::setPrecision(5);
	    printout(0, setw(3) << i);
	    printout(0, " " << setw(12) << real_norm);
	    printout(0, " " << setw(12) << imag_norm);
            TelePrompter::setPrecision(1);
	    printout(0, " " << setw(5) << rNodes);
	    printout(0, " " << setw(5) << iNodes);
	    printout(0, " " << setw(8) << dNorm_n);
	    printout(0, setw(9) << timer.getWallTime() << endl);	
	}
    }

    timer.stop();
    TelePrompter::printFooter(0, timer, 2);
    TelePrompter::setPrecision(oldprec);
}

void SCF::applyHelmholtzOperators_P(OrbitalVector &phi_np1,
                                  MatrixXd &F_n,
                                  OrbitalVector &phi_n,
                                  bool adjoint) {
    TelePrompter::printHeader(0, "Applying Helmholtz Operators");
    println(0, " Orb    OrbNorm       NormDiff       nNodes         Timing   ");
    TelePrompter::printSeparator(0, '-');
    int oldprec = TelePrompter::setPrecision(5);

    Timer timer;
    HelmholtzOperatorSet &H = *this->helmholtz;
    H.setPrecision(getOrbitalPrecision());

    phi_np1.clear();
    vector<Orbital *> arg_i_vec;

    int Ni = phi_n.size();
      OrbitalVector OrbVecChunk_i(0);//to store adresses of own i_orbs
      int OrbsIx[workOrbVecSize];//to store own orbital indices
      OrbitalVector rcvOrbs(0);//to store adresses of received orbitals
      int rcvOrbsIx[workOrbVecSize];//to store received orbital indices
      
      //make vector with adresses of own orbitals
      int i = 0;
      for (int Ix = MPI_rank; Ix < Ni; Ix += MPI_size) {
	OrbVecChunk_i.push_back(phi_n.getOrbital(Ix));//i orbitals
	OrbsIx[i++] = Ix;
      }

      Orbital *arg_i;
      bool first_iter = true;//part_1 is the Helmholtz part + the part2 accumulated in previous iterations
      for (int iter = 0;  iter >= 0; iter++) {
	//get a new chunk from other processes
	OrbVecChunk_i.getOrbVecChunk(OrbsIx, rcvOrbs, rcvOrbsIx, Ni, iter);
	Orbital *arg_i_1;
	for (int i = MPI_rank;  i<Ni ; i += MPI_size) {
	  Orbital &phi_i = phi_n.getOrbital(i);
	  if(first_iter){	    
	    arg_i_1 = getHelmholtzArgument_1(phi_i);
	  }else{
	    arg_i_1 = arg_i_vec[i/MPI_size];//set to part_1 + chunks so far
	  }
	  double coef_part1 = 1.0;
	  if(first_iter)coef_part1= -1.0/(2.0*pi);//only include factor once
	  arg_i = getHelmholtzArgument_2(i, rcvOrbsIx, F_n, rcvOrbs, arg_i_1, coef_part1, phi_i, adjoint);
	  if(first_iter){
	    arg_i_vec.push_back(arg_i);
	  }else{
	    arg_i_vec[i/MPI_size]=arg_i;
	  }
	}
	first_iter = false;
	rcvOrbs.clearVec(false);//reset to zero size orbital vector     
      }
      OrbVecChunk_i.clearVec(false);

      for (int i = MPI_rank;  i<Ni ; i += MPI_size) {
	  Timer timer;
	  Orbital &np1Phi_i = phi_np1.getOrbital(i);
	  arg_i = arg_i_vec[i/MPI_size];
	  H(i, np1Phi_i, *arg_i);
	  delete arg_i;
	
	int nNodes = np1Phi_i.getNNodes();

	Orbital &nPhi_i = phi_n.getOrbital(i);	  
	double norm_n = sqrt(nPhi_i.getSquareNorm());
	double norm_np1 = sqrt(np1Phi_i.getSquareNorm());
	double dNorm_n = fabs(norm_np1-norm_n);
	timer.stop();
	cout<< setw(3) << i;
	cout<< " " << setw(13) << norm_np1;
	cout<< " " << setw(13) << dNorm_n;
	cout<< " " << setw(8) << nNodes;
	cout<< setw(18) << timer.getWallTime() << endl;	
	//}
      }
 
    timer.stop();
    TelePrompter::printFooter(0, timer, 2);
    TelePrompter::setPrecision(oldprec);
}

Orbital* SCF::calcMatrixPart(int i, MatrixXd &M, OrbitalVector &phi) {
    vector<complex<double> > coefs;
    vector<Orbital *> orbs;

    int nOrbs = phi.size();
    for (int j = 0; j < nOrbs; j++) {
        double coef = M(i,j);
        // Linear scaling screening inserted here
        if (fabs(coef) > MachineZero) {
            Orbital &phi_j = phi.getOrbital(j);
            double norm_j = sqrt(phi_j.getSquareNorm());
            if (norm_j > 0.01*getOrbitalPrecision()) {
                coefs.push_back(coef);
                orbs.push_back(&phi_j);
            }
        }
    }

    Orbital &phi_i = phi.getOrbital(i);
    Orbital *result = new Orbital(phi_i);
    if (orbs.size() > 0) {
        Timer timer;
        this->add(*result, coefs, orbs, false);
        timer.stop();
        double time = timer.getWallTime();
        int nNodes = result->getNNodes();
        TelePrompter::printTree(2, "Added matrix part", nNodes, time);
    }
    return result;
}

Orbital* SCF::calcMatrixPart_P(int i_Orb, MatrixXd &M, OrbitalVector &phi) {
    vector<double> coefs;
    vector<Orbital *> orbs;
    Orbital &phi_i = phi.getOrbital(i_Orb);
    Orbital *result = new Orbital(phi_i);//should use workOrb?

    Timer timer;

    std::vector<complex<double> > U_Chunk;
    std::vector<Orbital *> orbChunk;
    int orbVecIx = 0;

    int nOrbs = phi.size();
    for (int iter = 0;  iter<MPI_size ; iter++) {
      int j_MPI=(MPI_size+iter-MPI_rank)%MPI_size;
      for (int j_Orb = j_MPI;  j_Orb < phi.size(); j_Orb+=MPI_size) {
	
	//    for (int j = 0; j < nOrbs; j++) {
	if(MPI_rank > j_MPI){
	  //send first bra, then receive ket
	  if (fabs(M(j_Orb,i_Orb)) > MachineZero)phi.getOrbital(i_Orb).send_Orbital(j_MPI, i_Orb);
	  if (fabs(M(i_Orb,j_Orb)) > MachineZero)workOrbVec.getOrbital(orbVecIx).Rcv_Orbital(j_MPI, j_Orb);
	}else if(MPI_rank < j_MPI){
	  if (fabs(M(i_Orb,j_Orb)) > MachineZero)workOrbVec.getOrbital(orbVecIx).Rcv_Orbital(j_MPI, j_Orb);
	  if (fabs(M(j_Orb,i_Orb)) > MachineZero)phi.getOrbital(i_Orb).send_Orbital(j_MPI, i_Orb);
	}else if(MPI_rank == j_MPI){//use phi directly
	}
	double coef = M(i_Orb,j_Orb);
	// Linear scaling screening inserted here
	if (fabs(coef) > MachineZero) {
	  //Orbital &phi_j = out.getOrbital(j_Orb);
	  double norm_j;
	  if(MPI_rank == j_MPI){
	    norm_j = sqrt(phi.getOrbital(j_Orb).getSquareNorm());
	  }else{
	    norm_j = sqrt(workOrbVec.getOrbital(orbVecIx).getSquareNorm());
	  }
	  if (norm_j > 0.01*getOrbitalPrecision()) {//NB: TODO should test before sending orbital!!
	    //coefs.push_back(coef);
	    //orbs.push_back(&workOrbVec.getOrbital(j_Orb));
	    
	    //push orbital in vector (chunk)
	    if(MPI_rank == j_MPI){
	      orbChunk.push_back(&phi.getOrbital(j_Orb));
	    }else{
	      orbChunk.push_back(&workOrbVec.getOrbital(orbVecIx++));
	    }
	    U_Chunk.push_back(coef);	
	  }
	  if(orbChunk.size()>=workOrbVecSize or iter >= MPI_size-1){
	    //Do the work for the chunk	  	  
	    this->add.inPlace(*result,U_Chunk, orbChunk, false);//can start with empty orbital
	    U_Chunk.clear();
	    orbChunk.clear();
	    orbVecIx = 0;
	  }
	}
      }
    }

    timer.stop();
    double time = timer.getWallTime();
    int nNodes = result->getNNodes();
    TelePrompter::printTree(2, "Added matrix part", nNodes, time);
    return result;
}
