
#include "qmfunctions.h"

namespace mrchem {

class OrbitalIterator final {
public:
    OrbitalIterator(OrbitalVector &Phi);

    bool next( bool symmetric = false , int max_recv = -1 );
    Orbital &get_orbital(int i) { return (this->received_orbitals)[i]; }
    int get_idx(int i) { return (this->received_orbital_index)[i]; }
    int get_size() { return this->received_orbitals.size(); }

protected:
    int iter;
    int received_counter; //number of orbitals fetched during this iteration
    int sent_counter; //number of orbitals sent during this iteration
    OrbitalVector received_orbitals;
    std::vector<int> received_orbital_index;
    OrbitalVector *orbitals;
    std::vector<std::vector<int>> orb_ix_map;//indices in the original orbital vector for each MPI
};

} //namespace mrchem
