#pragma once

#include <string>

#include "qmfunctions.h"

namespace mrchem {
class Molecule;

namespace gto_guess {
class OrbitalExp;

OrbitalVector initial_guess(double prec,
                            const Molecule &mol,
                            const std::string &bas_file,
                            const std::string &mo_file);
OrbitalVector initial_guess(double prec,
                            const Molecule &mol,
                            const std::string &bas_file,
                            const std::string &moa_file,
                            const std::string &mob_file);
void project(double prec, OrbitalVector &Phi, OrbitalExp &gto_exp);

} //namespace gto_guess
} //namespace mrchem
