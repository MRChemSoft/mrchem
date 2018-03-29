#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "ExchangePotential.h"
#include "Orbital.h"
#include "parallel.h"

using mrcpp::Printer;
using mrcpp::Timer;

//extern Orbital workOrb;
//extern OrbitalVector workOrbVec2;

namespace mrchem {

/** @brief constructor
 *
 * @param[in] P Poisson operator (does not take ownership)
 * @param[in] Phi vector of orbitals which define the exchange operator
 */
ExchangePotential::ExchangePotential(mrcpp::PoissonOperator &P, OrbitalVector &Phi, bool s)
            : screen(s),
              orbitals(&Phi),
              poisson(&P) {
    this->tot_norms = DoubleVector::Zero(Phi.size());
    this->part_norms = DoubleMatrix::Zero(Phi.size(), Phi.size());
}

/** @brief Prepare operator for application
 *
 * @param[in] prec reqested precision
 *
 * This will NOT precompute the internal exchange between the orbtials defining
 * the operator, which is done explicitly using setupInternal().
 */
void ExchangePotential::setup(double prec) {
    setApplyPrec(prec);
}

/** @brief Clears the Exchange Operator
 *
 *  Clears deletes the precomputed exchange contributions.
 */
void ExchangePotential::clear() {
    orbital::free(this->exchange);
    clearApplyPrec();
}

/** @brief Perform a unitary transormation among the precomputed exchange contributions
 *
 * @param[in] U unitary matrix defining the rotation
 */
void ExchangePotential::rotate(const ComplexMatrix &U) {
    if (this->exchange.size() == 0) return;

    OrbitalVector tmp = orbital::multiply(U, this->exchange, this->apply_prec);
    orbital::free(this->exchange);
    this->exchange = tmp;

    // NOTE: The following MPI point is currently NOT implemented!
    //
    // the last parameter, 1, means MPI will send only one orbital at a time
    // (because Exchange orbitals can be large for large molecules).
    // OrbitalAdder add(this->apply_prec, this->max_scale, 1);
    // add.rotate(this->exchange, U);
}

/** @brief Applies operator potential
 *
 *  @param[in] inp input orbital
 *
 * The exchange potential is applied to the given orbital. Checks first if this
 * particular exchange contribution has been precomputed.
 */
Orbital ExchangePotential::apply(Orbital inp) {
    if (this->apply_prec < 0.0) {
        MSG_ERROR("Uninitialized operator");
        return inp.paramCopy();
    }
    int i = testPreComputed(inp);
    if (i < 0) {
        println(1, "On-the-fly exchange");
        return calcExchange(inp);
    } else {
        println(1, "Precomputed exchange");
        return this->exchange[i].deepCopy();
    }
}

/** @brief Applies the adjoint of the operator
 *  \param[in] inp input orbital
 *
 * NOT IMPLEMENTED
 */
Orbital ExchangePotential::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

/** @brief Computes the exchange potential on the fly
 *
 *  \param[in] inp input orbital
 *
 * The exchange potential is computed and applied on the fly to the given orbital.
 */
Orbital ExchangePotential::calcExchange(Orbital phi_p) {
    Timer timer;

    double prec = this->apply_prec;
    OrbitalVector &Phi = *this->orbitals;
    mrcpp::PoissonOperator &P = *this->poisson;

    ComplexVector coef_vec(Phi.size());
    OrbitalVector orb_vec;

    for (int i = 0; i < Phi.size(); i++) {
        Orbital &phi_i = Phi[i];

        double spin_fac = getSpinFactor(phi_i.spin(), phi_p.spin());
        if (std::abs(spin_fac) < mrcpp::MachineZero) continue;

        // compute phi_ip = phi_i^dag * phi_p
        Orbital phi_ip = orbital::multiply(phi_i.dagger(), phi_p);

        // compute V_ip = P[phi_ip]
        Orbital V_ip = phi_p.paramCopy();
        if (phi_ip.hasReal()) {
            V_ip.alloc(NUMBER::Real);
            mrcpp::apply(prec, V_ip.real(), P, phi_ip.real());
        }
        if (phi_ip.hasImag()) {
            V_ip.alloc(NUMBER::Imag);
            mrcpp::apply(prec, V_ip.imag(), P, phi_ip.imag());
        }
        phi_ip.free();

        // compute phi_iip = phi_i * V_ip
        Orbital phi_iip = orbital::multiply(phi_i, V_ip);
        V_ip.free();

        coef_vec(i) = spin_fac/phi_i.squaredNorm();
        orb_vec.push_back(phi_iip);
        phi_iip.clear();
    }

    // compute ex_p = sum_i c_i*phi_iip
    Orbital ex_p = orbital::multiply(coef_vec, orb_vec, -1.0);
    orbital::free(orb_vec);

    timer.stop();
    double n = ex_p.getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(1, "Applied exchange", n, t);

    return ex_p;
/*
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
    */
}

/** @brief precomputes the exchange potential
 *
 *  @param[in] phi_p input orbital
 *
 * The exchange potential is (pre)computed among the orbitals that define the operator
 */
void ExchangePotential::setupInternal(double prec) {
    setApplyPrec(prec);
    if (mpi::orb_size > 1) NOT_IMPLEMENTED_ABORT;

    if (this->exchange.size() != 0) MSG_ERROR("Exchange not properly cleared");

    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &Ex = this->exchange;

    Timer timer;
    // Diagonal must come first because it's NOT in-place
    for (int i = 0; i < Phi.size(); i++) {
        calcInternal(i);
    }
    // Off-diagonal must come last because it IS in-place
    for (int i = 0; i < Phi.size(); i++) {
        for (int j = 0; j < i; j++) {
            calcInternal(i, j);
        }
    }

    int n = 0;
    // Collect info from the calculation
    for (int i = 0; i < Phi.size(); i++) {
        this->tot_norms(i) = Ex[i].norm();
        n += Ex[i].getNNodes();
    }

    timer.stop();
    double t = timer.getWallTime();
    Printer::printTree(0, "Hartree-Fock exchange", n, t);
/*
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
*/
}

/** @brief Computes the diagonal part of the internal exchange potential
 *
 *  \param[in] i orbital index
 *
 * The diagonal term K_ii is computed.
 */
void ExchangePotential::calcInternal(int i) {
    double prec = std::min(getScaledPrecision(i,i), 1.0e-1);

    Orbital &phi_i = (*this->orbitals)[i];
    mrcpp::PoissonOperator &P = *this->poisson;

    // compute phi_ii = phi_i^dag * phi_i
    Orbital phi_ii = orbital::multiply(phi_i.dagger(), phi_i, prec);

    // compute V_ii = P[phi_ii]
    Orbital V_ii = phi_i.paramCopy();
    if (phi_ii.hasReal()) {
        V_ii.alloc(NUMBER::Real);
        mrcpp::apply(prec, V_ii.real(), P, phi_ii.real());
    }
    if (phi_ii.hasImag()) {
        V_ii.alloc(NUMBER::Imag);
        mrcpp::apply(prec, V_ii.imag(), P, phi_ii.imag());
    }
    phi_ii.free();

    // compute phi_iii = phi_i * V_ii
    Orbital phi_iii = orbital::multiply(phi_i, V_ii, prec);
    phi_iii.rescale(1.0/phi_i.squaredNorm());
    this->part_norms(i,i) = phi_iii.norm();
    V_ii.free();

    this->exchange.push_back(phi_iii);
}

/** @brief computes the off-diagonal part of the exchange potential
 *
 *  \param[in] i first orbital index
 *  \param[in] j second orbital index
 *
 * The off-diagonal terms K_ij and K_ji are computed.
 */
void ExchangePotential::calcInternal(int i, int j) {
    mrcpp::PoissonOperator &P = *this->poisson;
    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &Ex = this->exchange;

    if (i == j) MSG_FATAL("Cannot handle diagonal term");
    if (Ex.size() != Phi.size()) MSG_FATAL("Size mismatch");
    if (Phi[i].hasImag() or Phi[j].hasImag()) MSG_FATAL("Orbitals must be real");

    double i_fac = getSpinFactor(Phi[i], Phi[j]);
    double j_fac = getSpinFactor(Phi[j], Phi[i]);

    double thrs = mrcpp::MachineZero;
    if (std::abs(i_fac) < thrs or std::abs(j_fac) < thrs) {
        this->part_norms(i,j) = 0.0;
        return;
    }

    // set correctly scaled precision for components ij and ji
    double prec = std::min(getScaledPrecision(i,j), getScaledPrecision(j,i));
    if (prec > 1.0e00) return;      // orbital does not contribute within the threshold
    prec = std::min(prec, 1.0e-1);  // very low precision does not work properly

    // compute phi_ij = phi_i^dag * phi_j (dagger NOT used, orbitals must be real!)
    Orbital phi_ij = orbital::multiply(Phi[i], Phi[j], prec);

    // compute V_ij = P[phi_ij]
    Orbital V_ij = Phi[i].paramCopy();
    if (phi_ij.hasReal()) {
        V_ij.alloc(NUMBER::Real);
        mrcpp::apply(prec, V_ij.real(), P, phi_ij.real());
    }
    if (phi_ij.hasImag()) {
        MSG_FATAL("Orbitals must be real");
        V_ij.alloc(NUMBER::Imag);
        mrcpp::apply(prec, V_ij.imag(), P, phi_ij.imag());
    }
    phi_ij.free();

    // compute phi_jij = phi_j * V_ij
    Orbital phi_jij = orbital::multiply(Phi[j], V_ij, prec);
    phi_jij.rescale(1.0/Phi[j].squaredNorm());
    this->part_norms(j,i) = phi_jij.norm();

    // compute phi_iij = phi_i * V_ij
    Orbital phi_iij = orbital::multiply(Phi[i], V_ij, prec);
    phi_jij.rescale(1.0/Phi[i].squaredNorm());
    this->part_norms(i,j) = phi_iij.norm();

    V_ij.free();

    // compute x_i += phi_jij
    Ex[i].add(i_fac, phi_jij);
    phi_jij.free();

    // compute x_j += phi_iij
    Ex[j].add(j_fac, phi_iij);
    phi_iij.free();
}

/** @brief computes the off-diagonal part of the exchange potential for own orbitals
 *  \param[in] i first orbital index
 *  \param[in] j second orbital index
 *
 * The off-diagonal term X_ij is computed.
 */
void ExchangePotential::calcInternal(int i, int j, Orbital &phi_i, Orbital &phi_j) {
    NOT_IMPLEMENTED_ABORT;
    /* waiting for MPI version
    //NB: this routine only compute ex_i, never ex_j
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
    */
}

/** @brief computes the off-diagonal part of the exchange potential for own orbitals
 *  and return V_ij for further use
 *  \param[in] i first orbital index
 *  \param[in] j second orbital index
 *
 * The off-diagonal term X_ij is computed.
 */
void ExchangePotential::calcInternal(int i, int j, Orbital &phi_i, Orbital &phi_j, Orbital* phi_iij) {
    NOT_IMPLEMENTED_ABORT;
    /* wait for MPI version
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
    */
}

/** @brief Test if a given contribution has been precomputed
 *
 * @param[in] phi_p orbital for which the check is performed
 *
 * If the given contribution has been precomputed, it is simply copied,
 * without additional recalculation.
 */
int ExchangePotential::testPreComputed(Orbital phi_p) const {
    const OrbitalVector &Phi = *this->orbitals;
    const OrbitalVector &Ex = this->exchange;

    int out = -1;
    if (Ex.size() == Phi.size()) {
        for (int i = 0; i < Phi.size(); i++) {
            if (&Phi[i].real() == &phi_p.real() and
                &Phi[i].imag() == &phi_p.imag()) {
                out = i;
                break;
            }
        }
    }
    return out;
}

/** @brief scale the relative precision based on norm
 *
 * The internal norms are saved between SCF iterations so that they can
 * be used to estimate the size of the different contributions to the total
 * exchange. The relative precision of the Poisson terms is scaled to give a
 * consistent _absolute_ pecision in the final output.
 */
double ExchangePotential::getScaledPrecision(int i, int j) const {
    double scaled_prec = this->apply_prec;
    if (this->screen) {
        double tNorm = this->tot_norms(i);
        double pNorm = std::max(this->part_norms(i,j), this->part_norms(j,i));
        if (tNorm > 0.0) scaled_prec *= tNorm/pNorm;
    }
    return scaled_prec;
}

/** @brief determines the exchange factor to be used in the calculation of the exact exchange
  *
  * @param [in] orb input orbital to which K is applied
  *
  * The factor is computed in terms of the occupancy of the two orbitals and in terms of the spin
  * 0.5 factors are used in order to preserve occupancy of the set of doubly occupied orbitals
  * this-> is the orbital defining the operator whereas the input orbital (orb) is the one
  * the operator is applied to
  *
  * Occupancy: Single/Double
  * Spin: alpha/beta
  *
  * K (this->) | orb (input) | factor
  * alpha      | alpha       | 1.0
  * alpha      | beta        | 0.0
  * alpha      | double      | 0.5
  * -------------------------------
  * beta       | alpha       | 0.0
  * beta       | beta        | 1.0
  * beta       | double      | 0.5
  * -------------------------------
  * double     | alpha       | 1.0
  * double     | beta        | 1.0
  * double     | double      | 1.0
  *
  */
double ExchangePotential::getSpinFactor(Orbital phi_i, Orbital phi_j) const {
    double out = 0.0;
         if (phi_j.spin() == SPIN::Paired) out = 1.0;
    else if (phi_i.spin() == SPIN::Paired) out = 0.5;
    else if (phi_i.spin() == phi_j.spin()) out = 1.0;
    return out;
}

} //namespace mrchem
