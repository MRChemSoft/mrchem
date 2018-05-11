#pragma once

#include "qmfunctions.h"

/** @file hydrogen_guess.h
 *
 * @brief Module for generating initial guess of hydrogen functions
 *
 * The hydrogen_guess namespace provides a single function in the public
 * interface, used to generate an initial guess of hydrogen eigenfunctions.
 * The initial guess requires no external input, but is not very accurate.
 */

namespace mrchem {
class Molecule;

namespace hydrogen_guess {
OrbitalVector initial_guess(double prec, const Molecule &mol, bool restricted, int zeta);
} //namespace hydrogen_guess

} //namespace mrchem
