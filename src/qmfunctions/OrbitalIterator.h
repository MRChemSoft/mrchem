
#include "qmfunctions.h"

namespace mrchem {

class OrbitalIterator final {
public:
    OrbitalIterator(OrbitalVector &Phi);

    bool next( bool symmetric = false );
    Orbital &get_orbital(int i) { return (this->received_orbitals)[i]; }
    int get_idx(int i) { return (this->received_orbital_index)[i]; }
    int get_size() { return this->received_orbitals.size(); }

protected:
    int iter;
    OrbitalVector received_orbitals;
    std::vector<int> received_orbital_index;
    OrbitalVector *orbitals;
};

} //namespace mrchem
