#pragma once

#include "qmfunctions/qmfunction_fwd.h"

/** @file sad.h
 *
 * @brief Module for generating initial guess as superposition of atomic densities
 *
 * The initial_guess::sad namespace provides functionality to setup an initial
 * guess from hydrogen functions and a superposition of atomic densities.
 */

namespace mrchem {
class Molecule;

namespace initial_guess {
namespace sad {

OrbitalVector setup(double prec, const Molecule &mol, bool restricted, int zeta);

} //namespace sad
} //namespace initial_guess
} //namespace mrchem
