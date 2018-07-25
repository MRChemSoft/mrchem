
#include "qmfunctions.h"

namespace mrchem {

class OrbitalIterator final {
public:
    OrbitalIterator(OrbitalVector &Phi);

    bool next();
    OrbitalChunk &get() { return this->chunk; }

protected:
    int iter;
    OrbitalChunk chunk;
    OrbitalVector *orbitals;
};

} //namespace mrchem
