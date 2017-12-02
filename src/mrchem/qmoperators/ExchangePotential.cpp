#include "ExchangePotential.h"
#include "OrbitalAdder.h"
#include "OrbitalMultiplier.h"
#include "PoissonOperator.h"
#include "MWConvolution.h"


using namespace std;
using namespace Eigen;

extern Orbital workOrb;
extern OrbitalVector workOrbVec2;

ExchangePotential::ExchangePotential(PoissonOperator &P,
                                     OrbitalVector &phi,
                                     double x_fac)
    : ExchangeOperator(P, phi, x_fac),
      exchange(phi) {
}

void ExchangePotential::setup(double prec) {
    setApplyPrec(prec);
    calcInternalExchange();
}

void ExchangePotential::clear() {
    this->exchange.clear();
    clearApplyPrec();
}

void ExchangePotential::rotate(MatrixXd &U) {
    // negative prec means precomputed exchange is cleared and should therefore not be rotated
    if (this->apply_prec > 0.0) {
        //the last parameter, 1, means MPI will send only one orbital at a time
	//(because Exchange orbitals can be large for large molecules). 
	OrbitalAdder add(this->apply_prec, this->max_scale, 1);
        add.rotate(this->exchange, U);
    }
}

Orbital* ExchangePotential::operator() (Orbital &phi_p) {
    Orbital *preOrb = testPreComputed(phi_p);
    if (preOrb != 0) {
        println(1, "Precomputed exchange");
        return preOrb;
    }
    return calcExchange(phi_p);
}

Orbital* ExchangePotential::adjoint(Orbital &phi_p) {
    NOT_IMPLEMENTED_ABORT;
}

/** Compute exchange on the fly */
Orbital* ExchangePotential::calcExchange(Orbital &phi_p) {
    Timer timer;

    OrbitalAdder add(this->apply_prec, this->max_scale);
    OrbitalMultiplier mult(this->apply_prec, this->max_scale);
    MWConvolution<3> apply(this->apply_prec, this->max_scale);
    PoissonOperator &P = *this->poisson;

    vector<complex<double> > coef_vec;
    vector<Orbital *> orb_vec;

    int nOrbs = this->orbitals->size();

    Orbital *ex_p = new Orbital(phi_p);

#ifdef HAVE_MPI

    OrbitalVector orbVecChunk_i(0); //to store adresses of own i_orbs
    OrbitalVector rcvOrbs(0);       //to store adresses of received orbitals
    vector<int> orbsIx;             //to store own orbital indices
    int rcvOrbsIx[workOrbVecSize];  //to store received orbital indices

    //make vector with adresses of own orbitals
    for (int Ix = mpiOrbRank;  Ix < nOrbs; Ix += mpiOrbSize) {
	orbVecChunk_i.push_back(this->orbitals->getOrbital(Ix));//i orbitals
	orbsIx.push_back(Ix);
    }

    for (int iter = 0;  iter >= 0; iter++) {
	//get a new chunk from other processes
	orbVecChunk_i.getOrbVecChunk(orbsIx, rcvOrbs, rcvOrbsIx, nOrbs, iter, workOrbVecSize, 2);
	for (int i = 0; i<rcvOrbs.size(); i++){
	    Orbital &phi_i = rcvOrbs.getOrbital(i);
	
	    double spinFactor = phi_i.getExchangeFactor(phi_p);
	    if (IS_EQUAL(spinFactor, 0.0)) continue;
	    Orbital *phi_ip = new Orbital(phi_p);
	    mult.adjoint(*phi_ip, 1.0, phi_i, phi_p);

	    Orbital *V_ip = new Orbital(phi_p);
	    if (phi_ip->hasReal()) {
		V_ip->allocReal();
		apply(V_ip->real(), P, phi_ip->real());
	    }
	    if (phi_ip->hasImag()) {
		V_ip->allocImag();
		apply(V_ip->imag(), P, phi_ip->imag());
	    }
	    delete phi_ip;

	    double multFac = - spinFactor * (this->x_factor / phi_i.getSquareNorm());
	    Orbital *phi_iip = new Orbital(phi_p);
	    mult(*phi_iip, multFac, phi_i, *V_ip);
	    delete V_ip;

	    coef_vec.push_back(1.0);
	    orb_vec.push_back(phi_iip);
	}
    }
    orbVecChunk_i.clearVec(false);
    rcvOrbs.clearVec(false);
    workOrbVec2.clear();
	
    add(*ex_p, coef_vec, orb_vec, true);

    for (int i = 0; i < orb_vec.size(); i++) delete orb_vec[i];
    orb_vec.clear();

#else
    for (int i = 0; i < nOrbs; i++) {
        Orbital &phi_i = this->orbitals->getOrbital(i);

        double spinFactor = phi_i.getExchangeFactor(phi_p);
        if (IS_EQUAL(spinFactor, 0.0)) continue;

        Orbital *phi_ip = new Orbital(phi_p);
        mult.adjoint(*phi_ip, 1.0, phi_i, phi_p);

        Orbital *V_ip = new Orbital(phi_p);
        if (phi_ip->hasReal()) {
            V_ip->allocReal();
            apply(V_ip->real(), P, phi_ip->real());
        }
        if (phi_ip->hasImag()) {
            V_ip->allocImag();
            apply(V_ip->imag(), P, phi_ip->imag());
        }
        delete phi_ip;

        double multFac = - spinFactor * (this->x_factor / phi_i.getSquareNorm());
        Orbital *phi_iip = new Orbital(phi_p);
        mult(*phi_iip, multFac, phi_i, *V_ip);
        delete V_ip;

        coef_vec.push_back(1.0);
        orb_vec.push_back(phi_iip);
    }

    add(*ex_p, coef_vec, orb_vec, true);

    for (int i = 0; i < orb_vec.size(); i++) delete orb_vec[i];
    orb_vec.clear();

#endif

    timer.stop();
    double n = ex_p->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Applied exchange", n, t);

    return ex_p;
}

/** Precompute the internal exchange */
void ExchangePotential::calcInternalExchange() {
    Timer timer;
    int nOrbs = this->orbitals->size();
    int n = 0;

    if (mpiOrbSize==1) {
	for (int i = 0; i < nOrbs; i++) {
	    calcInternal(i);
	    for (int j = 0; j < i; j++) {
		calcInternal(i,j);
	    }
	}
      
	for (int i = 0; i < nOrbs; i++) {
	    Orbital &ex_i = this->exchange.getOrbital(i);
	    this->tot_norms(i) = sqrt(ex_i.getSquareNorm());
	    n = max(n, ex_i.getNNodes());
	}
    } else {

#ifdef HAVE_MPI

	//use symmetri
	//send one orbital at a time
	MPI_Request request=MPI_REQUEST_NULL;
	MPI_Status status;
	//has to distribute the calculations evenly among processors
	OrbitalVector orbVecChunk_i(0); //to store adresses of own i_orbs
	OrbitalVector rcvOrbs(0);       //to store adresses of received orbitals
	vector<int> orbsIx;             //to store own orbital indices
	int rcvOrbsIx[workOrbVecSize];  //to store received orbital indices
      
	int sndtoMPI[workOrbVecSize];   //to store rank of MPI where orbitals were sent
	int sndOrbIx[workOrbVecSize];   //to store indices of where orbitals were sent
      
	//make vector with adresses of own orbitals
	for (int Ix = mpiOrbRank;  Ix < nOrbs; Ix += mpiOrbSize) {
	    orbVecChunk_i.push_back(this->orbitals->getOrbital(Ix));//i orbitals
	    orbsIx.push_back(Ix);
	}
      
	OrbitalAdder add(-1.0, this->max_scale);
	for (int i = 0; i < workOrbVecSize; i++) sndtoMPI[i]=-1;//init
	int mpiiter = 0;
	for (int iter = 0;  iter >= 0; iter++) {
	    mpiiter++;
	    //get a new chunk from other processes
	    sndOrbIx[0] = 0;//init
	    sndtoMPI[0] = -1;//init
	    orbVecChunk_i.getOrbVecChunk_sym(orbsIx, rcvOrbs, rcvOrbsIx, nOrbs, iter, sndtoMPI, sndOrbIx,1,2);
	    int rcv_left = 1;//normally we get one phi_jji per ii, maybe none
	    if (sndtoMPI[0] < 0 or sndtoMPI[0] == mpiOrbRank) rcv_left = 0;//we haven't sent anything, will not receive anything back	    
	    //convention: all indices with "i" are owned locally
	    for (int ii = 0; ii < orbVecChunk_i.size()+1 ; ii++){ //we may have to do one extra iteration to fetch all data
		int j = 0; //because we limited the size in getOrbVecChunk_sym to 1
		int i = ii;
		if (ii >= orbVecChunk_i.size()) i = 0; //so that phi_i is defined, but will not be used
	  
		Orbital &phi_i = orbVecChunk_i.getOrbital(i);
		Orbital *phi_iij = new Orbital(phi_i);
		int i_rcv = sndOrbIx[j];
		Orbital &phi_i_rcv = this->orbitals->getOrbital(i_rcv);
		Orbital *phi_jji_rcv = new Orbital(phi_i_rcv);	      
		if (rcvOrbs.size() > 0 and ii < orbVecChunk_i.size()) {
		    if (orbsIx[i] == rcvOrbsIx[j]) {
			//orbital should be own and i and j point to same orbital
			calcInternal(orbsIx[i]);
		    } else {
			Orbital &phi_j = rcvOrbs.getOrbital(j);	    
			if (rcvOrbsIx[j]%mpiOrbSize != mpiOrbRank) {
			    calcInternal(orbsIx[i], rcvOrbsIx[j], phi_i, phi_j, phi_iij);
			    //we send back the locally computed result to where j came from 
			    phi_iij->setOccupancy(orbsIx[i]);//We temporarily use Occupancy to send Orbital rank 
			    phi_iij->setSpin(orbVecChunk_i.size()-i-1);//We temporarily use Spin to send info about number of transfers left
			    phi_iij->setError( this->part_norms(rcvOrbsIx[j],orbsIx[i]));//We temporarily use Error to send part_norm
			    phi_iij->Isend_Orbital(rcvOrbsIx[j]%mpiOrbSize, mpiiter%10, request);
			} else {
			    //only compute j < i in own block 
			    if (rcvOrbsIx[j] < orbsIx[i]) calcInternal(orbsIx[i], rcvOrbsIx[j], phi_i, phi_j, phi_iij);
			}
		    }
		}
		if (rcv_left > 0) {
		    //we expect to receive data 
		    phi_jji_rcv->Rcv_Orbital(sndtoMPI[j], mpiiter%10);//we get back phi_jji from where we sent i
		}
		if (rcvOrbsIx[j]%mpiOrbSize != mpiOrbRank and rcvOrbs.size() > 0 and ii < orbVecChunk_i.size()){
		    MPI_Wait(&request, &status);//do not continue before isend is finished	      
		}
		if (rcv_left > 0) {
		    // phi_jji_rcv is ready
		    int i_rcv = sndOrbIx[j];
		    int j_rcv = phi_jji_rcv->getOccupancy();//occupancy was used to send rank!
		    rcv_left = phi_jji_rcv->getSpin();//occupancy was used to send number of transfers left!
		    if (i_rcv!=j_rcv) {
			phi_jji_rcv->setOccupancy(this->orbitals->getOrbital(j_rcv).getOccupancy());//restablish occupancy		
			phi_jji_rcv->setSpin(this->orbitals->getOrbital(j_rcv).getSpin());//restablish spin		
			this->part_norms(i_rcv,j_rcv) = phi_jji_rcv->getError();//does not seem to matter?
			phi_jji_rcv->setError(this->orbitals->getOrbital(j_rcv).getError());//restablish error		
			//j_rcv and i_rcv are now the active indices j_rcv is the remote orbital
			// compute x_i_rcv += phi_jji_rcv j_rcv is the remote orbital here i_rcv is what has been used to compute phi_jji_rcv
			this->part_norms(j_rcv,i_rcv) = sqrt(phi_jji_rcv->getSquareNorm());
			Orbital &ex_i_rcv = this->exchange.getOrbital(i_rcv);
			Orbital &phi_j_rcv = this->orbitals->getOrbital(j_rcv);
			double i_factor_rcv = phi_i_rcv.getExchangeFactor(phi_j_rcv);
			add.inPlace(ex_i_rcv, i_factor_rcv, *phi_jji_rcv);
		    }
		}
	  
		if (phi_jji_rcv != 0) delete phi_jji_rcv;
		if (phi_iij != 0) delete phi_iij;
	    }
	    rcvOrbs.clearVec(false);//reset to zero size orbital vector     	
	}
	orbVecChunk_i.clearVec(false);
	workOrbVec2.clear();

	for (int i = 0; i < nOrbs; i++) {
	    if (i%mpiOrbSize == mpiOrbRank) {
		Orbital &ex_i = this->exchange.getOrbital(i);
		this->tot_norms(i) = sqrt(ex_i.getSquareNorm());
		n = max(n, ex_i.getNNodes());
	    } else {
		this->tot_norms(i) = 0.0;
	    }
	}
	//tot_norms are used for screening. Since we use symmetri, we might need factors from others
	MPI_Allreduce(MPI_IN_PLACE, &this->tot_norms(0), nOrbs, MPI_DOUBLE, MPI_SUM, mpiCommOrb);

#endif
    }
   
    timer.stop();
    double t = timer.getWallTime();
    TelePrompter::printTree(0, "Hartree-Fock exchange", n, t);
}

void ExchangePotential::calcInternal(int i) {
    OrbitalAdder add(-1.0, this->max_scale);
    OrbitalMultiplier mult(-1.0, this->max_scale);
    MWConvolution<3> apply(-1.0, this->max_scale);

    PoissonOperator &P = *this->poisson;
    Orbital &phi_i = this->orbitals->getOrbital(i);

    double prec = getScaledPrecision(i, i);
    prec = min(prec, 1.0e-1);

    mult.setPrecision(prec);
    apply.setPrecision(prec);

    // compute phi_ii = phi_i^dag * phi_i
    Orbital *phi_ii = new Orbital(phi_i);
    mult.adjoint(*phi_ii, 1.0, phi_i, phi_i);

    // compute V_ii = P[phi_ii]
    Orbital *V_ii = new Orbital(phi_i);
    if (phi_ii->hasReal()) {
        V_ii->allocReal();
        apply(V_ii->real(), P, phi_ii->real());
    }
    if (phi_ii->hasImag()) {
        V_ii->allocImag();
        apply(V_ii->imag(), P, phi_ii->imag());
    }
    if (phi_ii != 0) delete phi_ii;

    double fac_iii = -(this->x_factor/phi_i.getSquareNorm());

    // compute phi_iii = phi_i * V_ii
    Orbital *phi_iii = new Orbital(phi_i);
    mult(*phi_iii, fac_iii, phi_i, *V_ii);
    this->part_norms(i,i) = sqrt(phi_iii->getSquareNorm());
    if (V_ii != 0) delete V_ii;

    // compute x_i += phi_iii
    Orbital &ex_i = this->exchange.getOrbital(i);
    add.inPlace(ex_i, 1.0, *phi_iii);
    if (phi_iii != 0) delete phi_iii;
}

void ExchangePotential::calcInternal(int i, int j) {
    OrbitalAdder add(-1.0, this->max_scale);
    OrbitalMultiplier mult(-1.0, this->max_scale);
    MWConvolution<3> apply(-1.0, this->max_scale);

    PoissonOperator &P = *this->poisson;
    Orbital &phi_i = this->orbitals->getOrbital(i);
    Orbital &phi_j = this->orbitals->getOrbital(j);

    double i_factor = phi_i.getExchangeFactor(phi_j);
    double j_factor = phi_j.getExchangeFactor(phi_i);
    if (IS_EQUAL(i_factor, 0.0) or IS_EQUAL(j_factor, 0.0)) {
        this->part_norms(i,j) = 0.0;
        return;
    }

    double prec = std::min(getScaledPrecision(i, j), getScaledPrecision(j, i));

    if (prec > 1.0e00) return;
    prec = min(prec, 1.0e-1);

    mult.setPrecision(prec);
    apply.setPrecision(prec);

    // compute phi_ij = phi_i^dag * phi_j
    Orbital *phi_ij = new Orbital(phi_i);
    mult.adjoint(*phi_ij, 1.0, phi_i, phi_j);

    // compute V_ij = P[phi_ij]
    Orbital *V_ij = new Orbital(phi_i);
    if (phi_ij->hasReal()) {
        V_ij->allocReal();
        apply(V_ij->real(), P, phi_ij->real());
    }
    if (phi_ij->hasImag()) {
        V_ij->allocImag();
        apply(V_ij->imag(), P, phi_ij->imag());
    }
    if (phi_ij != 0) delete phi_ij;

    // compute phi_jij = phi_j * V_ij
    double fac_jij = -(this->x_factor/phi_j.getSquareNorm());
    Orbital *phi_jij = new Orbital(phi_i);
    mult(*phi_jij, fac_jij, *V_ij, phi_j);
    this->part_norms(j,i) = sqrt(phi_jij->getSquareNorm());

    // compute phi_iij = phi_i * V_ij
    double fac_iij = -(this->x_factor/phi_i.getSquareNorm());
    Orbital *phi_iij = new Orbital(phi_i);
    mult.adjoint(*phi_iij, fac_iij, *V_ij, phi_i);
    this->part_norms(i,j) = sqrt(phi_iij->getSquareNorm());

    if (V_ij != 0) delete V_ij;

    // compute x_i += phi_jij
    Orbital &ex_i = this->exchange.getOrbital(i);
    add.inPlace(ex_i, i_factor, *phi_jij);
    if (phi_jij != 0) delete phi_jij;

    // compute x_j += phi_iij
    Orbital &ex_j = this->exchange.getOrbital(j);
    add.inPlace(ex_j, j_factor, *phi_iij);
    if (phi_iij != 0) delete phi_iij;
}

//NB: this routine only compute ex_i, never ex_j
void ExchangePotential::calcInternal(int i, int j, Orbital &phi_i, Orbital &phi_j) {
    OrbitalAdder add(-1.0, this->max_scale);
    OrbitalMultiplier mult(-1.0, this->max_scale);
    MWConvolution<3> apply(-1.0, this->max_scale);

    PoissonOperator &P = *this->poisson;

    double i_factor = phi_i.getExchangeFactor(phi_j);
    double j_factor = phi_j.getExchangeFactor(phi_i);
    if (IS_EQUAL(i_factor, 0.0) or IS_EQUAL(j_factor, 0.0)) {
        this->part_norms(i,j) = 0.0;
        return;
    }

    double prec = std::min(getScaledPrecision(i, j), getScaledPrecision(j, i));

    if (prec > 1.0e00) return;
    prec = min(prec, 1.0e-1);

    mult.setPrecision(prec);
    apply.setPrecision(prec);

    // compute phi_ij = phi_i^dag * phi_j
    Orbital *phi_ij = new Orbital(phi_i);
    mult.adjoint(*phi_ij, 1.0, phi_i, phi_j);

    // compute V_ij = P[phi_ij]
    Orbital *V_ij = new Orbital(phi_i);
    if (phi_ij->hasReal()) {
        V_ij->allocReal();
        apply(V_ij->real(), P, phi_ij->real());
    }
    if (phi_ij->hasImag()) {
        V_ij->allocImag();
        apply(V_ij->imag(), P, phi_ij->imag());
    }
    if (phi_ij != 0) delete phi_ij;

    // compute phi_jij = phi_j * V_ij
    double fac_jij = -(this->x_factor/phi_j.getSquareNorm());
    Orbital *phi_jij = new Orbital(phi_i);
    mult(*phi_jij, fac_jij, *V_ij, phi_j);
    this->part_norms(j,i) = sqrt(phi_jij->getSquareNorm());

    //part_norms(i,j) MUST be computed to use for symmetric screening
    // compute phi_iij = phi_i * V_ij
    Orbital *phi_iij = new Orbital(phi_i);
    double fac_iij = -(this->x_factor/phi_i.getSquareNorm());
    mult.adjoint(*phi_iij, fac_iij, *V_ij, phi_i);
    this->part_norms(i,j) = sqrt(phi_iij->getSquareNorm());
    if (phi_iij != 0) delete phi_iij;

    if (V_ij != 0) delete V_ij;

    // compute x_i += phi_jij
    Orbital &ex_i = this->exchange.getOrbital(i);
    add.inPlace(ex_i, i_factor, *phi_jij);
    if (phi_jij != 0) delete phi_jij;
}

void ExchangePotential::calcInternal(int i, int j, Orbital &phi_i, Orbital &phi_j, Orbital* phi_iij) {
    //to check: adjoint or not?
    OrbitalAdder add(-1.0, this->max_scale);
    OrbitalMultiplier mult(-1.0, this->max_scale);
    MWConvolution<3> apply(-1.0, this->max_scale);

    PoissonOperator &P = *this->poisson;

    double i_factor = phi_i.getExchangeFactor(phi_j);
    double j_factor = phi_j.getExchangeFactor(phi_i);
    if (IS_EQUAL(i_factor, 0.0) or IS_EQUAL(j_factor, 0.0)) {
        this->part_norms(i,j) = 0.0;
        return;
    }

    double prec = std::min(getScaledPrecision(i, j), getScaledPrecision(j, i));
    if (prec > 1.0e00) return;
    prec = min(prec, 1.0e-1);

    mult.setPrecision(prec);
    apply.setPrecision(prec);

    // compute phi_ij = phi_i^dag * phi_j
    Orbital *phi_ij = new Orbital(phi_i);
    mult.adjoint(*phi_ij, 1.0, phi_i, phi_j);

    // compute V_ij = P[phi_ij]
    Orbital *V_ij = new Orbital(phi_i);
    if (phi_ij->hasReal()) {
        V_ij->allocReal();
        apply(V_ij->real(), P, phi_ij->real());
    }
    if (phi_ij->hasImag()) {
        V_ij->allocImag();
        apply(V_ij->imag(), P, phi_ij->imag());
    }
    if (phi_ij != 0) delete phi_ij;

    // compute phi_jij = phi_j * V_ij
    double fac_jij = -(this->x_factor/phi_j.getSquareNorm());
    Orbital *phi_jij = new Orbital(phi_i);
    mult(*phi_jij, fac_jij, *V_ij, phi_j);
    this->part_norms(j,i) = sqrt(phi_jij->getSquareNorm());

    // compute x_i += phi_jij
    Orbital &ex_i = this->exchange.getOrbital(i);
    add.inPlace(ex_i, i_factor, *phi_jij);
    if (phi_jij != 0) delete phi_jij;

    //part_norms(i,j) MUST be computed to use for symmetric screening
    // compute phi_iij = phi_i * V_ij
    double fac_iij = -(this->x_factor/phi_i.getSquareNorm());
    mult.adjoint(*phi_iij, fac_iij, *V_ij, phi_i);
    this->part_norms(i,j) = sqrt(phi_iij->getSquareNorm());

    if (V_ij != 0) delete V_ij;

    //This part (phi_iij) is sent to owner of j orbital if it is on another MPI
    if (j%mpiOrbSize == mpiOrbRank) {
	// compute x_j += phi_iij
        Orbital &ex_j = this->exchange.getOrbital(j);
        add.inPlace(ex_j, j_factor, *phi_iij);
    }
}

Orbital* ExchangePotential::testPreComputed(Orbital &phi_p) {
    int nOrbs = this->orbitals->size();
    
    for (int i = 0; i < nOrbs; i++) {
        Orbital &phi_i  = this->orbitals->getOrbital(i);
        Orbital *ex_i = this->exchange[i];
        if (&phi_i.real() == &phi_p.real() and &phi_i.imag() == &phi_p.imag() and ex_i->getNNodes() != 0) {
            Orbital *result = new Orbital(phi_p);
            // Deep copy of orbital
            OrbitalAdder add(this->apply_prec, this->max_scale);
            add.inPlace(*result, 1.0, *ex_i);
            return result;
        }
    }
    return 0;
}
