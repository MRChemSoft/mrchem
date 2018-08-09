#include "MRCPP/Printer"

#include "OrbitalIterator.h"
#include "Orbital.h"
#include "parallel.h"

namespace mrchem {

OrbitalIterator::OrbitalIterator(OrbitalVector &Phi)
        : iter(0),
          orbitals(&Phi) {
}

// Receive all the orbitals, for non symmetric treatment
// Receive half of the orbitals, for symmetric treatment
// Send own orbitals to half of the MPI, for symmetric treatment
// The algorithm ensures that at each iteration, each MPI is communicating
// with only one other MPI. This ensures good communication parallelism.
bool OrbitalIterator::next( bool symmetric ) {
    mpi::free_foreign(*this->orbitals); //delete foreign orbitals
    //NB: for now we completely delete orbitals at each iteration. Could reuse the memory in next iteration.
    this->received_orbital_index.clear();
    this->received_orbitals.clear();

    if (this->iter >= mpi::orb_size) {
        // We have talked to all MPIs -> return
        this->iter = 0;
        return false;
    }

    //  We communicate with one MPI per iteration:
    //   - Symmetric: do two iterations, send in one step receive in the other
    //   - send ALL my orbitals to other_rank
    //   - receive ALL orbitals owned by other_rank
    //   - other_rank increases by one for each iteration, modulo size
    // Each MPI pair will be treated exactly once at the end of the iteration process

    int my_rank = mpi::orb_rank;
    int max_rank = mpi::orb_size;

    int nstep = 1;// send and receive in the same iteration
    if (symmetric) nstep = 2;// one iteration is receive, one iteration is send

    for (int step = 0; step < nstep; step++) {

	if(this->iter >= mpi::orb_size) break;

	// Figure out which MPI to talk to
	int other_rank = (max_rank + this->iter - my_rank)%max_rank;
	//note: my_rank= (max_rank + this->iter - other_rank)%max_rank

	if (other_rank == my_rank) {
	    // We send/receive to ourself
	    for (int i = 0; i < this->orbitals->size(); i++) {
		Orbital &phi_i = (*this->orbitals)[i];
		if (mpi::my_orb(phi_i)) {
		    this->received_orbital_index.push_back(i);
		    this->received_orbitals.push_back(phi_i);
		}
	    }
	} else {

	    bool IamReceiver = true;//default, both send and receive
	    bool IamSender   = true;//default, both send and receive

	    if(symmetric) {
		// Either send or receive
		// Figure out which one will send and receive
		//(my_rank + other_rank )%2 flips each time other_rank changes by one
		//(my_rank < other_rank) ensures that my_rank and other_rank are opposite
		if ( (my_rank + other_rank + (my_rank < other_rank))%2 ) {
		    IamReceiver = false;
		} else {
		    IamSender   = false;
		}
	    }

	    for (int i = 0; i < this->orbitals->size(); i++) {
		int tag = max_rank*this->iter + i;
		Orbital &phi_i = (*this->orbitals)[i];
		int phi_rank = phi_i.rankID();
		if (phi_rank == other_rank and IamReceiver) {
		    mpi::recv_orbital(phi_i, other_rank, tag);
		    this->received_orbital_index.push_back(i);
		    this->received_orbitals.push_back(phi_i);
		}
		if (phi_rank == my_rank and IamSender) {
		    mpi::send_orbital(phi_i, other_rank, tag);
		}
	    }
	}
	this->iter++;
    }
    return true;
}

} //namespace mrchem

