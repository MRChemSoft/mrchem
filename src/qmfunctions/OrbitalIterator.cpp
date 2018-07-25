#include "MRCPP/Printer"

#include "OrbitalIterator.h"
#include "Orbital.h"
#include "parallel.h"

namespace mrchem {

OrbitalIterator::OrbitalIterator(OrbitalVector &Phi)
        : iter(0),
          orbitals(&Phi) {
}

bool OrbitalIterator::next() {
    mpi::free_foreign(*this->orbitals);
    this->chunk.clear();

    if (this->iter >= mpi::orb_size) {
        // We have talked to all MPIs -> return
        this->iter = 0;
        return false;
    }


    //  We talk to one MPI at a time:
    //   - send ALL my orbitals to other_rank
    //   - receive ALL orbitals owned by other_rank
    //   - other_rank increases by one for each iteration, modulo size
    //  All orbitals belonging to one MPI will be sent exactly
    //  once per iteration of next()

    int my_rank = mpi::orb_rank;
    int max_rank = mpi::orb_size;

    // Figure out which MPI to talk to
    int other_rank = (max_rank + this->iter - my_rank)%max_rank;
    //note: my_rank = (max_rank + this->iter - other_rank)%max_rank

    if (other_rank  == my_rank) {
	// We talk to ourselves
	for (int i = 0; i < this->orbitals->size(); i++) {
	    Orbital &phi_i = (*this->orbitals)[i];
	    if (mpi::my_orb(phi_i)) {
		this->chunk.push_back(std::make_tuple(i, phi_i));
	    }
	}
    } else {
	for (int i = 0; i < this->orbitals->size(); i++) {
	    int tag = max_rank*this->iter + i;
	    Orbital &phi_i = (*this->orbitals)[i];
	    int phi_rank = phi_i.rankID();
	    if (phi_rank == other_rank) {
		mpi::recv_orbital(phi_i, other_rank, tag);
		this->chunk.push_back(std::make_tuple(i, phi_i));
	    }
	    if (phi_rank == my_rank) {
		mpi::send_orbital(phi_i, other_rank, tag);
	    }
	}
    }

    this->iter++;
    return true;
}

} //namespace mrchem

