#include "MRCPP/Printer"

#include "parallel.h"

#include "OrbitalVector.h"

namespace mrchem {

/** OrbitalVector constructor
 *
 * ne is number of electrons.
 * mult is spin multiplicity.
 * rest is spin restricted.
 *
 * All orbital functions are uninitialized.
 *
 */
OrbitalVector::OrbitalVector(int ne, int mult, bool rest)
        : in_use(false) {
    int de = ne - (mult - 1);
    if (de%2 != 0)  MSG_ERROR("Invalid multiplicity");

    int na, nb, nd;
    if (rest) {
        na = mult - 1;
        nb = 0;
        nd = de/2;
    } else {
        na = de/2 + (mult - 1);
        nb = de/2;
        nd = 0;
    }

    for (int i = 0; i < nd; i++) push_back(SPIN::Paired);
    for (int i = 0; i < na; i++) push_back(SPIN::Alpha);
    for (int i = 0; i < nb; i++) push_back(SPIN::Beta);
}

/** Copy constructor
 *
 * New orbitals are constructed as shallow copies of the input set.
 *
 */
OrbitalVector::OrbitalVector(const OrbitalVector &inp_vec)
        : in_use(false) {
    for (int i = 0; i < inp_vec.size(); i++) {
        this->orbitals.push_back(inp_vec[i]);
    }
}

/** Assignment operator
 *
 * New orbitals are constructed as shallow copies of the input set.
 *
 */
OrbitalVector& OrbitalVector::operator=(const OrbitalVector &inp_vec) {
    this->clear();
    if (this != &inp_vec) {
        for (int i = 0; i < inp_vec.size(); i++) {
            Orbital out_i;
            out_i = inp_vec[i];
            this->push_back(out_i);
        }
    }
    return *this;
}

/** Deep copy
 *
 * New orbitals are constructed as deep copies of the input set.
 *
 */
OrbitalVector OrbitalVector::deepCopy() {
    OrbitalVector &inp_vec = *this;
    OrbitalVector out_set;
    for (int i = 0; i < inp_vec.size(); i++) {
        Orbital out_i = inp_vec[i].deepCopy();
        out_set.push_back(out_i);
    }
    return out_set;
}

/** Parameter copy
 *
 * New orbitals are constructed as parameter copies of the input set.
 *
 */
OrbitalVector OrbitalVector::paramCopy() const {
    const OrbitalVector &inp_vec = *this;
    OrbitalVector out_set;
    for (int i = 0; i < inp_vec.size(); i++) {
        Orbital out_i = inp_vec[i].paramCopy();
        out_set.push_back(out_i);
    }
    return out_set;
}

/** Clears the orbital vector
 *
 * Leaves an empty vector. The orbitals are not freed.
 *
 */
void OrbitalVector::clear() {
    this->orbitals.clear();
}

/** Frees each orbital in the vector
 *
 * Leaves an empty vector. Orbitals are freed.
 *
 */
void OrbitalVector::free() {
    for (int i = 0; i < this->size(); i++) {
        Orbital &orb = (*this)[i];
        orb.free();
    }
    this->orbitals.clear();
}

/** Adjoin two vectors
 *
 * The orbitals of the input vector are appended to
 * (*this) vector, the ownership is transferred. Leaves
 * the input vector empty.
 *
 */
void OrbitalVector::adjoin(OrbitalVector &inp) {
    for (int i = 0; i < inp.size(); i++) {
        this->push_back(inp[i]);
    }
    inp.clear();
}

/** Disjoin vector in two parts
 *
 * All orbitals of a particular spin is collected in a new vector
 * and returned. These orbitals are removed from (*this) vector,
 * and the ownership is transferred.
 *
 */
OrbitalVector OrbitalVector::disjoin(int spin) {
    OrbitalVector out;
    OrbitalVector tmp;
    for (int i = 0; i < this->size(); i++) {
        Orbital &orb_i = (*this)[i];
        if (orb_i.spin() == spin) {
            out.push_back(orb_i);
        } else {
            tmp.push_back(orb_i);
        }
    }
    this->clear();
    *this = tmp;
    return out;
}


/** Returns the number of occupied orbitals */
int OrbitalVector::getNOccupied() const {
    int nOcc = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital &orb = (*this)[i];
        if (orb.occ() > 0) nOcc++;
    }
    return nOcc;
}

/** Returns the number of empty orbitals */
int OrbitalVector::getNEmpty() const {
    int nEmpty = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital &orb = (*this)[i];
        if (orb.occ() == 0) nEmpty++;
    }
    return nEmpty;
}

/** Returns the number of singly occupied orbitals */
int OrbitalVector::getNSingly() const {
    int nSingly = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital &orb = (*this)[i];
        if (orb.occ() == 1) nSingly++;
    }
    return nSingly;
}

/** Returns the number of doubly occupied orbitals */
int OrbitalVector::getNDoubly() const {
    int nDoubly = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital &orb = (*this)[i];
        if (orb.occ() == 1) nDoubly++;
    }
    return nDoubly;
}

/** Returns the number of paired orbitals */
int OrbitalVector::getNPaired() const {
    int nPaired = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital &orb = (*this)[i];
        if (orb.spin() == SPIN::Paired) nPaired++;
    }
    return nPaired;
}

/** Returns the number of alpha orbitals */
int OrbitalVector::getNAlpha() const {
    int nAlpha = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital &orb = (*this)[i];
        if (orb.spin() == SPIN::Alpha) nAlpha++;
    }
    return nAlpha;
}

/** Returns the number of beta orbitals */
int OrbitalVector::getNBeta() const {
    int nBeta = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital &orb = (*this)[i];
        if (orb.spin() == SPIN::Beta) nBeta++;
    }
    return nBeta;
}

/** Returns the number of electrons with the given spin
 *
 * Paired spin (default input) returns the total number of electrons.
 *
 */
int OrbitalVector::getNElectrons(int spin) const {
    int nElectrons = 0;
    for (int i = 0; i < this->size(); i++) {
        const Orbital &orb = (*this)[i];
        if (spin == SPIN::Paired) {
            nElectrons += orb.occ();
        } else if (spin == SPIN::Alpha) {
            if (orb.spin() == SPIN::Paired or orb.spin() == SPIN::Alpha) {
                nElectrons += 1;
            }
        } else if (spin == SPIN::Beta) {
            if (orb.spin() == SPIN::Paired or orb.spin() == SPIN::Beta) {
                nElectrons += 1;
            }
        } else {
            MSG_ERROR("Invalid spin argument");
        }
    }
    return nElectrons;
}

/** Returns the spin multiplicity of the vector */
int OrbitalVector::getMultiplicity() const {
    int nAlpha = getNElectrons(SPIN::Alpha);
    int nBeta = getNElectrons(SPIN::Beta);
    int S = abs(nAlpha - nBeta);
    return S + 1;
}

/** Returns a vector containing the orbital errors */
DoubleVector OrbitalVector::getErrors() const {
    int nOrbs = this->size();
    DoubleVector errors = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        const Orbital &orb = (*this)[i];
        errors(i) = orb.error();
    }
    return errors;
}

/** Assign errors to each orbital.
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void OrbitalVector::setErrors(const DoubleVector &errors) {
    int nOrbs = this->size();
    if (nOrbs != errors.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < nOrbs; i++) {
        Orbital &orb = (*this)[i];
        orb.setError(errors(i));
    }
}

/** Returns a vector containing the orbital spins */
IntVector OrbitalVector::getSpins() const {
    int nOrbs = this->size();
    IntVector spins = IntVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        const Orbital &orb = (*this)[i];
        spins(i) = orb.spin();
    }
    return spins;
}

/** Assigns spin to each orbital
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void OrbitalVector::setSpins(const IntVector &spins) {
    int nOrbs = this->size();
    if (nOrbs != spins.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < nOrbs; i++) {
        Orbital &orb = (*this)[i];
        orb.setSpin(spins(i));
    }
}

/** Returns a vector containing the orbital occupancies */
IntVector OrbitalVector::getOccupancies() const {
    int nOrbs = this->size();
    IntVector occ = IntVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        const Orbital &orb = (*this)[i];
        occ(i) = orb.occ();
    }
    return occ;
}

/** Assigns spin to each orbital
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void OrbitalVector::setOccupancies(const IntVector &occ) {
    int nOrbs = this->size();
    if (nOrbs != occ.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < nOrbs; i++) {
        Orbital &orb = (*this)[i];
        orb.setOcc(occ(i));
    }
}

/** Returns a vector containing the orbital square norms */
DoubleVector OrbitalVector::getSquaredNorms() const {
    int nOrbs = this->size();
    DoubleVector norms = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        const Orbital &orb = (*this)[i];
        norms(i) = orb.squaredNorm();
    }
    return norms;
}

/** Returns a vector containing the orbital norms */
DoubleVector OrbitalVector::getNorms() const {
    int nOrbs = this->size();
    DoubleVector norms = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        const Orbital &orb = (*this)[i];
        norms(i) = orb.norm();
    }
    return norms;
}

/** Normalize all orbitals in the set */
void OrbitalVector::normalize() {
    for (int i = 0; i < this->size(); i++) {
        Orbital &orb = (*this)[i];
        orb.normalize();
    }
}

/** Gram-Schmidt orthogonalize orbitals within the set */
void OrbitalVector::orthogonalize() {
    for (int i = 0; i < this->size(); i++) {
        Orbital &orb_i = (*this)[i];
        for (int j = 0; j < i; j++) {
            Orbital &orb_j = (*this)[j];
            orb_i.orthogonalize(orb_j);
        }
    }
}

/** Orthogonalize the out orbital against all orbitals in inp */
void OrbitalVector::orthogonalize(OrbitalVector inp_vec) {
    for (int i = 0; i < this->size(); i++) {
        Orbital &orb_i = (*this)[i];
        orb_i.orthogonalize(inp_vec);
    }
}

/** In place rotation (unitary transformation) of orbitals */
void OrbitalVector::rotate(const ComplexMatrix &U, double prec) {
    if (mpi::orb_size > 1) NOT_IMPLEMENTED_ABORT;

    OrbitalVector out;
    for (int i = 0; i < U.rows(); i++) {
        const ComplexVector &c = U.row(i);
	Orbital out_i = orbital::add(c, *this, prec);
        out.push_back(out_i);
    }
    this->free();
    *this = out;
}

/** Collective printing of all orbitals in the vector */
std::ostream& OrbitalVector::print(std::ostream &o) const {
    int oldPrec = mrcpp::Printer::setPrecision(15);
    mrcpp::Printer::setScientific();
    o << "*OrbitalVector: ";
    o << std::setw(4) << this->size()          << " orbitals  ";
    o << std::setw(4) << this->getNOccupied()  << " occupied  ";
    o << std::setw(4) << this->getNElectrons() << " electrons " << std::endl;
    o << "------------------------------";
    o << "------------------------------\n";
    o << "   n    sqNorm               Occ Spin  Error\n";
    o << "------------------------------";
    o << "------------------------------\n";
    for (int i = 0; i < this->size(); i++) {
        const Orbital &orb = (*this)[i];
        o << std::setw(4) << i << orb;
    }
    o << "------------------------------";
    o << "------------------------------\n";
    mrcpp::Printer::setPrecision(oldPrec);
    return o;
}

/*
struct Metadata{
    int Norbitals;
    int spin[workOrbVecSize];
    int occupancy[workOrbVecSize];
    int NchunksReal[workOrbVecSize];
    int NchunksImag[workOrbVecSize];
    int Ix[workOrbVecSize];
    double error[workOrbVecSize];
};
*/
//send an orbitalvector with MPI
/*
void OrbitalVector::send_OrbVec(int dest, int tag, vector<int> &orbsIx, int start, int maxcount){
#ifdef HAVE_MPI
    MPI_Status status;
    MPI_Comm comm=mpiCommOrb;

    Metadata Orbinfo;

    Orbinfo.Norbitals = min(this->size() - start, maxcount);

    Orbital* orb_i;
    for (int i = start; i < this->size() && (i-start<maxcount); i++) {
	int i_out=i-start;
	orb_i = &this->getOrbital(i);
	Orbinfo.spin[i_out] = orb_i->getSpin();
	Orbinfo.occupancy[i_out] = orb_i->getOccupancy();
	Orbinfo.error[i_out] = orb_i->getError();
	if(orb_i->hasReal()){
	    Orbinfo.NchunksReal[i_out] = orb_i->real().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
	}else{Orbinfo.NchunksReal[i_out] = 0;}
	if(orb_i->hasImag()){
	    Orbinfo.NchunksImag[i_out] = orb_i->imag().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
	}else{Orbinfo.NchunksImag[i_out] = 0;}
	Orbinfo.Ix[i_out] = orbsIx[i];
    }

    int count=sizeof(Metadata);
    MPI_Send(&Orbinfo, count, MPI_BYTE, dest, tag, comm);

    for (int i = start; i <  this->size() && (i-start<maxcount); i++) {
	int i_out=i-start;
	orb_i = &this->getOrbital(i);
	if(orb_i->hasReal())Send_SerialTree(&orb_i->real(), Orbinfo.NchunksReal[i_out], dest, 2*i_out+1+tag, comm);
	if(orb_i->hasImag())Send_SerialTree(&orb_i->imag(), Orbinfo.NchunksImag[i_out], dest, 2*i_out+2+tag, comm);
    }

#endif
}
*/

/*
#ifdef HAVE_MPI
//non-blocking send an orbitalvector with MPI
void OrbitalVector::Isend_OrbVec(int dest, int tag, vector<int> &orbsIx, int start, int maxcount, MPI_Request &request){
    MPI_Status status;
    MPI_Comm comm=mpiCommOrb;

    Metadata Orbinfo;

    Orbinfo.Norbitals = min(this->size() - start, maxcount);

    Orbital* orb_i;
    for (int i = start; i < this->size() && (i-start<maxcount); i++) {
	int i_out=i-start;
	orb_i = &this->getOrbital(i);
	Orbinfo.spin[i_out] = orb_i->getSpin();
	Orbinfo.occupancy[i_out] = orb_i->getOccupancy();
	Orbinfo.error[i_out] = orb_i->getError();
	if(orb_i->hasReal()){
	    Orbinfo.NchunksReal[i_out] = orb_i->real().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
	}else{Orbinfo.NchunksReal[i_out] = 0;}
	if(orb_i->hasImag()){
	    Orbinfo.NchunksImag[i_out] = orb_i->imag().getSerialFunctionTree()->nodeChunks.size();//should reduce to actual number of chunks
	}else{Orbinfo.NchunksImag[i_out] = 0;}
	Orbinfo.Ix[i_out] = orbsIx[i];
    }

    int count=sizeof(Metadata);

    MPI_Isend(&Orbinfo, count, MPI_BYTE, dest, tag, comm, &request);

    for (int i = start; i <  this->size() && (i-start<maxcount); i++) {
	int i_out=i-start;
	orb_i = &this->getOrbital(i);
	if(orb_i->hasReal())ISend_SerialTree(&orb_i->real(), Orbinfo.NchunksReal[i_out], dest, 2*i_out+1+tag, comm, request);
	if(orb_i->hasImag())ISend_SerialTree(&orb_i->imag(), Orbinfo.NchunksImag[i_out], dest, 2*i_out+2+tag, comm, request);
    }
}
#endif
*/

/*
//receive an orbitalvector with MPI
void OrbitalVector::Rcv_OrbVec(int source, int tag, int *orbsIx, int& workOrbVecIx){
#ifdef HAVE_MPI
    MPI_Status status;
    MPI_Comm comm=mpiCommOrb;
    OrbitalVector* workOrbVec_p = 0;//pointer to the workOrbVec to use

    if(workOrbVec2.inUse){
	workOrbVec_p = &workOrbVec2;//use workOrbVec2
    }else{
	if(workOrbVec.inUse){
	    workOrbVec_p = &workOrbVec;//use workOrbVec2 (workOrbVec may be in use already)
	}else{
	    MSG_ERROR("Error did not find which WorkOrbVec to use");
	}
    }


    Metadata Orbinfo;

    int count=sizeof(Metadata);
    MPI_Recv(&Orbinfo, count, MPI_BYTE, source, tag, comm, &status);

    Orbital* orb_i;
    for (int i = 0; i < Orbinfo.Norbitals; i++) {
	if(this->size() < workOrbVecIx+1){
	    //make place for the incoming orbital adresses
	    Orbital *orb = new Orbital(Orbinfo.occupancy[i], Orbinfo.spin[i]);
	    this->orbitals.push_back(orb);
	}
	//the orbitals are stored in workOrbVec, and only the metadata/adresses is copied into OrbitalVector
	orb_i = &workOrbVec_p->getOrbital(workOrbVecIx);

	orb_i->setSpin(Orbinfo.spin[i]);
	orb_i->setOccupancy(Orbinfo.occupancy[i]);
	orb_i->setError(Orbinfo.error[i]);

	if(Orbinfo.NchunksReal[i]>0){
	    if(not orb_i->hasReal()){
		//We must have a tree defined for receiving nodes. Define one:
		orb_i->allocReal();
	    }
	    Rcv_SerialTree(&orb_i->real(), Orbinfo.NchunksReal[i], source, 2*i+1+tag, comm);}

	if(Orbinfo.NchunksImag[i]>0){
	    if(not orb_i->hasImag()){
		//We must have a tree defined for receiving nodes. Define one:
		orb_i->allocImag();
	    }
	    Rcv_SerialTree(&orb_i->imag(), Orbinfo.NchunksImag[i], source, 2*i+2+tag, comm);
	}else{
	    //&(this->imag())=0;
	}
	orbsIx[workOrbVecIx] = Orbinfo.Ix[i];
	//the orbitals are stored in workOrbVec, and only the metadata/adresses is copied into OrbitalVector
	this->getOrbital(workOrbVecIx).shallowCopy(workOrbVec_p->getOrbital(workOrbVecIx));
	workOrbVecIx++;
    }

#endif
}
*/

/** Send and receive a chunk of an OrbitalVector
 * this : Vector with orbitals locally (owned by this MPI process)
 * myOrbsIx : indices of myOrbs in the orbital vector
 * rcvOrbs : Vector with orbitals received by from other processes
 * rcvOrbsIx : indices of rcvOrbs in the orbital vector
 * size : total number of Orbitals in OrbitalVector (all MPI. this->size is for only one MPI)
 * iter0 : iteration to start with, to know which chunk is to be treated
 * NB: it is assumed that the orbitals are evenly distributed among MPI processes (or as evenly as possible)
 */
/*
void OrbitalVector::getOrbVecChunk(vector<int> &myOrbsIx, OrbitalVector &rcvOrbs, int* rcvOrbsIx, int size, int& iter0, int maxOrbs_in, int workIx){

    int maxOrbs = workOrbVecSize;//max number of orbital to send or receive per iteration
    if(maxOrbs_in>0)maxOrbs = maxOrbs_in;

    int maxsizeperOrbvec = (size + mpiOrbSize-1)/mpiOrbSize;

    int RcvOrbVecIx = 0;
*/
    /* many index scales:
     * Orbital vector index. max = size
     * MPI_iter index. The iteration count through all MPI processes. max = mpiOrbSize
     * Index of local orbital (owned by the local MPI process). max = maxsizeperOrbvec
     * Chunk iteration. max = (maxsizeperOrbvec*mpiOrbSize +  maxOrbs-1)/maxOrbs
     * For accounting, we assume that all MPI are filled with maxsizeperOrbvec;
     * this is in order to know where to restart
     * Last iteration is for cleanups: clear send and receive vectors and set workOrbVec flag as available.
     */
/*
    int MPI_iter0 = (iter0*maxOrbs)/maxsizeperOrbvec;//which MPI_iter to start with
    int start = (iter0*maxOrbs)%maxsizeperOrbvec;//where to start in the first iteration
    int start0 = start;//save

    rcvOrbs.clearVec(false);//we restart from beginning of vector

    workOrbVec.inUse = true;
    if(workIx!=0){
	workOrbVec2.inUse = true;
    }else{
	workOrbVec2.inUse = false;
    }

    if(MPI_iter0>=mpiOrbSize){
	iter0=-2;//to indicate that we are finished with receiving AND sending all orbitals
    }

    for (int MPI_iter = MPI_iter0;  true; MPI_iter++) {
	int maxcount = maxOrbs;//place left first iteration
	if (MPI_iter > MPI_iter0) maxcount = maxOrbs+start0-(MPI_iter-MPI_iter0)*maxsizeperOrbvec;//place left after first iteration
	if(maxcount<=0) break;
	if(MPI_iter>=mpiOrbSize){
	    iter0=-2;//to indicate that we are finished with receiving AND sending all orbitals
	    break;
	}
	int rcv_MPI = (mpiOrbSize+MPI_iter-mpiOrbRank)%mpiOrbSize;//rank of process to communicate with
	if(mpiOrbRank > rcv_MPI){
	    //send first, then receive
	    this->send_OrbVec(rcv_MPI, MPI_iter, myOrbsIx, start, maxcount);
	    rcvOrbs.Rcv_OrbVec(rcv_MPI, MPI_iter, rcvOrbsIx, RcvOrbVecIx);
	}else if(mpiOrbRank < rcv_MPI){
	    //receive first, then send
	    rcvOrbs.Rcv_OrbVec(rcv_MPI,MPI_iter, rcvOrbsIx, RcvOrbVecIx);
	    this->send_OrbVec(rcv_MPI, MPI_iter, myOrbsIx, start, maxcount);
	}else{
	    for (int i = start;  i < this->size() && (i-start<maxcount) ; i++){
		rcvOrbsIx[RcvOrbVecIx] = myOrbsIx[i];
		if(rcvOrbs.size()<=RcvOrbVecIx){
		    rcvOrbs.push_back(this->getOrbital(i));
		}else{
		    rcvOrbs.getOrbital(RcvOrbVecIx) = this->getOrbital(i);
		}
		RcvOrbVecIx++;
	    }
	}
	start = 0;//always start at 0 after first iteration
    }
}
*/

/** Send and receive a chunk of an OrbitalVector
 * Assumes that a symmetric operator is calculated, so that only
 * half of the orbitals need to be sent. i.e exactly one of (i,j)or(j,i) is sent.
 * The orbitals are sent and received from different processors.
 * this : Vector with orbitals locally (owned by this MPI process)
 * myOrbsIx : indices of myOrbs in the orbital vector (in)
 * rcvOrbs : Vector with orbitals received by from other processes (out)
 * rcvOrbsIx : indices of rcvOrbs in the orbital vector (out)
 * size : total number of Orbitals in OrbitalVector (in)
 * iter0 : iteration to start with, to know which chunk is to be treated (inout)
 * sndtoMPI : rank of MPI where the orbitals have been sent (optional out)
 * sndOrbIx : index of orbitals that have been sent (optional out)
 * workIx : set if need to use another workOrbVec (optional in)
 * NB: it is assumed that the orbitals are evenly distributed among MPI processes (or as evenly as possible)
 */
/*
void OrbitalVector::getOrbVecChunk_sym(vector<int> &myOrbsIx, OrbitalVector &rcvOrbs, int* rcvOrbsIx, int size, int& iter0, int* sndtoMPI, int* sndOrbIx, int maxOrbs_in, int workIx){

    int maxOrbs = workOrbVecSize;//max number of orbital to send or receive per iteration
    if(maxOrbs_in>0)maxOrbs = maxOrbs_in;

    int maxsizeperOrbvec = (size + mpiOrbSize-1)/mpiOrbSize;

    int RcvOrbVecIx = 0;
#ifdef HAVE_MPI

    MPI_Request request[mpiOrbSize];
    for (int i = 0; i < mpiOrbSize; i++) request[i] = MPI_REQUEST_NULL;
    int nRequests = 0;
    MPI_Status status[mpiOrbSize];
*/
    /* many index scales:
     * Orbital vector index. max = size
     * MPI_iter index. The iteration count through all MPI processes. max = mpiOrbSize
     * Index of local orbital (owned by the local MPI process). max = maxsizeperOrbvec
     * Chunk iteration. max = (maxsizeperOrbvec*mpiOrbSize + maxOrbs-1)/maxOrbs
     * For accounting, we assume that all MPI are filled with maxsizeperOrbvec;
     * this is in order to know where to restart
     */
/*
    int MPI_iter0 = (iter0*maxOrbs)/maxsizeperOrbvec;//which MPI_iter to start with
    int start = (iter0*maxOrbs)%maxsizeperOrbvec;//where to start in the first iteration
    int start0 = start;//save
    rcvOrbs.clearVec(false);//we restart from beginning of vector

    workOrbVec.inUse = true;
    if(workIx!=0){
	workOrbVec2.inUse = true;
    }else{
	workOrbVec2.inUse = false;
    }

    if(MPI_iter0 >= (mpiOrbSize/2 + 1) ){
	iter0=-2;//to indicate that we are finished with receiving AND sending all orbitals
    }

    int isnd=0;
    for (int MPI_iter = MPI_iter0;  MPI_iter < mpiOrbSize; MPI_iter++) {
	int maxcount = maxOrbs;//place left
	if (MPI_iter > MPI_iter0) maxcount = maxOrbs+start0-(MPI_iter-MPI_iter0)*maxsizeperOrbvec;//place left
	if(maxcount <= 0)  break;//chunk is full
	if(MPI_iter >= (mpiOrbSize/2 + 1) ){
	    iter0=-2;//to indicate that we are finished with receiving AND sending all orbitals
	    break;
	}
	int rcv_MPI = abs((mpiOrbSize-MPI_iter+mpiOrbRank)%mpiOrbSize);//rank of process to receive from
	int snd_MPI = (mpiOrbSize+MPI_iter+mpiOrbRank)%mpiOrbSize;//rank of process to send to

	if(rcv_MPI == snd_MPI){
	    //only one of them do need to do the job: the one with lowest rank (of course!)
	    if(mpiOrbRank>rcv_MPI)rcv_MPI = -1;//to indicate receive nothing
	    if(mpiOrbRank<rcv_MPI)snd_MPI = -1;//to indicate send nothing
	}
	if(mpiOrbRank != rcv_MPI){
	    //non-blocking send first, then receive
	    if(snd_MPI>=0){
		for (int i = start;  i < this->size() && (i-start<maxcount) && sndOrbIx!=0 && sndOrbIx!=0; i++){
		    sndtoMPI[isnd]=snd_MPI;
		    sndOrbIx[isnd]=myOrbsIx[i];
		    isnd++;
		}
		this->Isend_OrbVec(snd_MPI, MPI_iter, myOrbsIx, start, maxcount, request[nRequests++]);
	    }
	    if(rcv_MPI>=0){
		rcvOrbs.Rcv_OrbVec(rcv_MPI, MPI_iter, rcvOrbsIx, RcvOrbVecIx);
	    }
	}else{
	    //own orbitals
	    for (int i = start;  i < this->size() && (i-start<maxcount) ; i++){
		rcvOrbsIx[RcvOrbVecIx] = myOrbsIx[i];
		if(rcvOrbs.size()<=RcvOrbVecIx){
		    rcvOrbs.push_back(this->getOrbital(i));
		}else{
		    rcvOrbs.getOrbital(RcvOrbVecIx) = this->getOrbital(i);
		}
		RcvOrbVecIx++;
	    }
	}
	start = 0;//always start at 0 after first iteration
    }
    MPI_Waitall(nRequests,request, status);
#endif
}
*/

} //namespace mrchem
