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

    if (this->iter == 0) {
        // First we talk to ourselves
        for (int i = 0; i < this->orbitals->size(); i++) {
            Orbital &phi_i = (*this->orbitals)[i];
            if (mpi::my_orb(phi_i)) {
                this->chunk.push_back(std::make_tuple(i, phi_i));
            }
        } 
    } else {
        // Then we talk to two MPIs at the time (can be the same):
        //   - send ALL my orbitals to send_rank
        //   - receive ALL orbitals owned by recv_rank
        //   - send_rank decreases by one for each iteration of next()
        //   - recv_rank increases by one for each iteration of next()
        int my_rank = mpi::orb_rank;
        int max_rank = mpi::orb_size;

        // Figure out which MPI to talk to
        int send_rank = my_rank - this->iter;
        int recv_rank = my_rank + this->iter;

        // Wrap around round robin
        if (send_rank < 0)         send_rank += max_rank;
        if (recv_rank >= max_rank) recv_rank -= max_rank;

        // Each orbital will be sent exactly once per iteration of next()
        for (int i = 0; i < this->orbitals->size(); i++) {
            int tag = 100000*this->iter + i;
            Orbital &phi_i = (*this->orbitals)[i];
            int phi_rank = phi_i.rankID();
            if (phi_rank == recv_rank) {
                mpi::recv_orbital(phi_i, recv_rank, tag);
                this->chunk.push_back(std::make_tuple(i, phi_i));
            }
            if (phi_rank == my_rank) {
                mpi::send_orbital(phi_i, send_rank, tag);
            }
        }
    }
    this->iter++;
    return true;
}

} //namespace mrchem

