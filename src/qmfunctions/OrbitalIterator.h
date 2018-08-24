
#include "qmfunctions.h"

namespace mrchem {

class OrbitalIterator final {
public:
    OrbitalIterator(OrbitalVector &Phi);

    bool next( bool symmetric = false , int max_recv = -1 );
    int get_size() { return this->received_orbitals.size(); }
    Orbital &get_orbital(int i) { return (this->received_orbitals)[i]; }
    int get_idx(int i) { return (this->received_orbital_index)[i]; }
    int get_idx_sent(int i) { return (this->sent_orbital_index)[i]; }
    int get_mpirank_sent(int i) { return (this->sent_orbital_mpirank)[i]; }
    int get_step(int i) { return (this->rcv_step)[i]; }
    int get_sent_size() { return this->sent_orbital_index.size(); }
    int get_iter() { return this->iter; }

protected:
    int iter;
    int received_counter; //number of orbitals fetched during this iteration
    int sent_counter; //number of orbitals sent during this iteration
    OrbitalVector received_orbitals;
    std::vector<int> received_orbital_index; //index of the orbitals received
    std::vector<int> sent_orbital_index; //indices of the orbitals sent out
    std::vector<int> sent_orbital_mpirank; //mpi rank (=rankID) of the orbitals sent out
    std::vector<int> rcv_step; //for symmetric treatment: if the send was at step=0 or 1
    OrbitalVector *orbitals;
    std::vector<std::vector<int>> orb_ix_map;//indices in the original orbital vector for each MPI
};

} //namespace mrchem
