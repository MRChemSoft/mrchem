#include "qmfunctions.h"

namespace mrchem {
class Molecule;
class Nuclei;

namespace hydrogen_guess {
extern int PT[29][2];

OrbitalVector initial_guess(double prec, const Molecule &mol, bool restricted, int zeta);
OrbitalVector project(double prec, const Nuclei &nucs, int zeta);
void populate(OrbitalVector &vec, int N, int spin);

} //namespace hydrogen_guess
} //namespace mrchem
