#pragma once

#include "mrchem.h"

class Molecule;
class Nuclei;

namespace mrchem {

double compute_nuclear_repulsion(const Nuclei &nucs);

} //namespace mrchem
