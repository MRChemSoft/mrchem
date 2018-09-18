#pragma once

#include "qmfunctions/qmfunction_fwd.h"

/** @file core.h
 *
 * @brief Module for generating initial guess of hydrogen functions
 *
 * The initial_guess::core namespace provides functionality to setup an
 * initial guess of hydrogen eigenfunctions.
 */

namespace mrchem {
class Molecule;
class Nuclei;

namespace initial_guess {
namespace core {

OrbitalVector setup(double prec, const Molecule &mol, bool restricted, int zeta);
OrbitalVector project_ao(double prec, const Nuclei &nucs, int spin, int zeta);

} //namespace core
} //namespace initial_guess
} //namespace mrchem
