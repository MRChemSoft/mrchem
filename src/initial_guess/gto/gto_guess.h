#pragma once

#include <string>

#include "qmfunctions.h"

/** @file gto_guess.h
 *
 * @brief Module for generating initial guess of GTO molecular orbitals
 *
 * The gto_guess namespace provides a single function in the public interface,
 * used to generate an initial guess of GTO molecular orbitals. The initial guess
 * requires external input to provide basis set information and MO coefficients.
 */

namespace mrchem {
class Molecule;

namespace gto_guess {
OrbitalVector initial_guess(double prec,
                            const Molecule &mol,
                            const std::string &bas_file,
                            const std::string &mo_file);
OrbitalVector initial_guess(double prec,
                            const Molecule &mol,
                            const std::string &bas_file,
                            const std::string &moa_file,
                            const std::string &mob_file);
} //namespace gto_guess
} //namespace mrchem
