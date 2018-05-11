#pragma once

#include "qmfunctions.h"

/** @file core_guess.h
 *
 * @brief Module for generating initial guess of hydrogen functions
 *
 * The core_guess namespace provides a single function in the public
 * interface, used to generate an initial guess of hydrogen eigenfunctions.
 * The initial guess requires no external input, but is not very accurate.
 */

namespace mrchem {
class Molecule;

namespace core_guess {
OrbitalVector initial_guess(double prec, const Molecule &mol, bool restricted, int zeta);
} //namespace core_guess

} //namespace mrchem
